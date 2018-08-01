#!/usr/bin/perl -s

#	perl -w -s genome_maker.pl
#	-sample_pileup_file=sample.mpileup
#	-reference_mapped_consensus=sample.vs.reference.consensus
#	-lastz_analysed_file=sample.lastz_analysed_file
#       -logfile=sample.genome_maker.log
#	-ref_correct_start=790
#	-ref_correct_stop=9417
#	-contigs=sample.fasta
#	-output=sample.genome.fasta


## Very old notes.  Instead of inserting coding genome start & stops
## we now simply trim UTR's off the reference sequences - saves Perl 
## some work!

#  NB  We'll just use hxb2.newtopline.fasta for the reference_mapped_consensus option for now.
#  Can input the reference-based derived corrected consensus in there later on if bwa/stampy have problems.
#  ref_correct_start and ref_correct_stop refer to the Gag start and Nef stop positions (i.e. coding regions).
#  These should be loaded from a separate file perhaps - the unused variable "reference_name" could be used
#  to look them up from a table.  We have at most 200 or so reference names and coords so will take ~0seconds
#  to load this in at run-time :)



open(LOG,">$logfile");

## load in contigs so can work out contig length, from which we can work out depths of coverage for reverse complement alignments!?!?
open(CONTIGS,"$contigs") || die $!;
while(<CONTIGS>){
	chomp($_);
	if(/>(\S+)/){
		$contigname = $1;
	}
	else{
		$contigseqs{$contigname} .= $_;
	}
}
close(CONTIGS) || die $!;

foreach $contig(keys %contigseqs){
	$contiglengths{$contig} = length($contigseqs{$contig});
}

### Firstly load in mpileup so have depths at each position from each contig
open(SAMPLE_PILEUP,"$sample_pileup_file") || die $!;
while(<SAMPLE_PILEUP>){
	if(/^(\S+)\s+(\S+)\s+(\S+)\s+(\S+)/){
		$contig = $1;
		$pos = $2;
		$base = $3;
		$depth = $4;
		$depth_hash{$contig}{$pos} = $depth;
		$revcompcontig = $contig."_RC";
		$revcomppos = $contiglengths{$contig} - $pos +1;
		$depth_hash{$revcompcontig}{$revcomppos} = $depth;
	}
}
close(SAMPLE_PILEUP) || die $!;

### Secondly load in reference-based consensus (from mapping reads to closest reference and correcting the consensus by
### majority (or highest frequency e.g. 40%A 30%G 30%C would give an A) variants.
### Will use this sequence to "fill in the gaps" where IVA assembly (or other e.g. SPADES) has failed.
### NB Currently simply inserting uncorrected reference sequence as it saves a round of bwa/stampy and (hopefully) any
### majority errors can be corrected in the genome-based bwa/stampy consensus correction step further down the pipeline.

open(REF_SEQ,"$reference_mapped_consensus") || die $!;
while(<REF_SEQ>){
	if(/>(\S+)/){
# UNUSED	$reference_name = $1;
	}
	else{
		chomp($_);
		$reference_mapped_consensus_sequence.= uc($_);
	}
}
close(REF_SEQ) || die $!;

for($i = 1; $i <=length($reference_mapped_consensus_sequence) ; $i++){
	$reference_mapped_consensus_map{$i} = substr($reference_mapped_consensus_sequence,($i-1),1);
}
###could just use: $pos = 0;while($reference_mapped_consensus_sequence =~ /(\S)/g){$pos++;$map{$pos} = $1;}

# take in output of lastz_analyser.pl - technically could run this within lastz_analyser.pl script - kept separate to check intermediary files and for validation right now...
# only want the part with the reference<->contig mapping - output now adjusted to be just this - all the alignment and indel/mutation info goes into log.lastz_analyse logfile
open(LASTZ_ANALYSED,"$lastz_analysed_file") || die $!;
$oldblurb = "";
while(<LASTZ_ANALYSED>){
	if(/^(\d+)\t(\S)\t(\S+)(.*)\n/){
		$referencepos = $1;
		$referencebase = $2;
		$sample = $3;
		$blurb = $4;
		if($referencepos == 1){
			$samplepos = 0;
		}
		if($blurb ne ""){
			if(/PROBLEM_OVERLAP/){
				# use the contig with the greatest depth at that position...
				%tempseqs = ();
				while($blurb =~ /con(\S+)\s+Length:(\d+)\s+AlignQuality:(\d+\.*\d*)\s+(\S+)\s+(\d+)\t(\S+)\s*/g){
					$contignumber = $1;
					$length = $2;
					$alignquality = $3;
					$fullcontigname = $4;
					$pos = $5;
					$base = $6;
					$depth = $depth_hash{$fullcontigname}{$pos} // 0;
					#  This is a tricky one - two contigs overlapping = easy problem.  
					#  But three contigs overlapping: 
					#  con1: 10000 bases G 
					#  con2: 4000 bases C 
					#  con3: 7000 bases C
					#  Should we go with the maximum (on the off chance contigs 2 and 3 are iffy assemblies, likely to be true if Con1 a lot longer than Con2 and Con3)
					#  Or sum up so C=11000 bases > 10000 G?
					#  Also get some instances where the same contig has 2 alignments - aaagh!  Here we should probably use alignment quality score - this should sort out
					#  secondary alignments (e.g. first alignment is ~9000nt long, second is 300nt long for a particular region in env.  Throw away second alignment...
					$tempseqs{$fullcontigname}{$alignquality}{$depth}{$base}++;
				}
				%problembase = ();
				foreach $contig (keys %tempseqs){
						@tempaqs = reverse sort {$a <=> $b} (keys %{$tempseqs{$contig}});
						$maxaq = $tempaqs[0];
						print LOG "Contig:$contig\tmaxaq:$maxaq\t";
						@tempdepths = reverse sort {$a <=> $b} (keys %{$tempseqs{$contig}{$maxaq}});
						$maxdepth = $tempdepths[0];
						print LOG "maxdepth:$maxdepth\n";
						# this picks the final base for a contig at random if the alignquality and depth are all the same
						# the occurrence of a base in @contigfinalbase is determined by the number of alignments
						# if aq and depth same have no way of knowing which is better
						# but ups the chances of that base on the off-chance (!) there's a particularly horrendous set of
						# alignments all with the same aq and depth...
					foreach $base (keys %{$tempseqs{$contig}{$maxaq}{$maxdepth}}){
						$problembase{$base} += $maxdepth;
					}
				} # end foreach $contig keys %tempseqs
				$fullmaxdepth = 0;
				foreach $base (keys %problembase){
					if($problembase{$base} > $fullmaxdepth){
						$fullmaxdepth = $problembase{$base};
					}
				}
				@problembase = ();
				foreach $base (keys %problembase){
					if($problembase{$base} == $fullmaxdepth){
						push(@problembase,$base);
					}
				}
				if(scalar @problembase == 1){
					$problembase = $problembase[0];
					$problemdepth = $fullmaxdepth;
				}
				else{
					$i = int rand(scalar @problembase);
					$problembase = $problembase[$i];
					$problemdepth = $fullmaxdepth;
				}
				$sampleseq{$sample} .= $problembase;
				print LOG "Final base:$problembase\tFinal depth:$problemdepth\n";
				$samplepos++;
				$samplemap{$sample}{$samplepos} = $referencepos;
			}
			elsif(/MUTATION:\S\d+(\S)/){
				$sampleseq{$sample} .= $1;
				$samplepos++;
				$samplemap{$sample}{$samplepos} = $referencepos;
			}
			elsif(/DELETION_IN_ALIGNMENT_TO_REFERENCE/){
			# $sampleseq{$sample} .= '-';
			# COMMENTED OUT TO SAVE NEEDING A SED TO REMOVE HYPHENS FROM SEQUENCE LATER
			# WAS USEFUL WHEN CHECKING CONTIGS WERE BEING CORRECTLY JOINED IN VALIDATION STEP...
			}
			elsif(/NO_ALIGNMENT_TO_REFERENCE_FOR_THIS_POSITION/){
				$samplepos++;
				$samplemap{$sample}{$samplepos} = $referencepos;
				if(($referencepos>=$ref_correct_start) && ($referencepos <= $ref_correct_stop)){
					print LOG "$reference_mapped_consensus inserted $reference_mapped_consensus_map{$referencepos} at refpos: $referencepos for $sample at samplepos: $samplepos\n";
					$sampleseq{$sample} .= $reference_mapped_consensus_map{$referencepos};
				}
				else{
					print LOG "N inserted as outside CDS boundaries (".$ref_correct_start." - ".$ref_correct_stop.") at refpos: $referencepos for $sample at: $samplepos\n";
					$sampleseq{$sample} .= 'N';
				}
			}
			elsif($blurb ne ""){
				$samplepos++;
				$samplemap{$sample}{$samplepos} = $referencepos;
			$sampleseq{$sample} .= $referencebase;
			}
		$oldblurb = $blurb;
		} # end if blurb ne ""
		else{
			#### shouldn't need samplepos++ here - maybe add a Z instead of - below to see if anything falls through the net...
			$sampleseq{$sample} .= "-";
		}
	} # end main reference <-> sample map
	elsif(/INS_\d+\t(\S+)\t(.*)/){
		$sample = $1;
		$blurb = $2;
		%no_insertion_depth = ();
		$max_no_insertion_depth = 0;
		%probleminsertion = ();
		# this line checks whether there are contigs which DON'T have an insertion - we should probably use the depth of these as a guide to whether we really
		# have an insertion or not - i.e. don't start putting in tons of insertions based on one or more dodgy contigs...
		while($oldblurb =~ /con(\S+)\s+Length:(\d+)\s+AlignQuality:(\d+\.*\d*)\s+(\S+)\s+(\d+)\t(\S+)\s*/g){
			$contignumber = $1;
			$length = $2;
			$alignquality = $3;
			$fullcontigname = $4;
			$pos = $5;
			$base = $6;
			$depth = $depth_hash{$fullcontigname}{$pos} // 0;
			$no_insertion_depth{$fullcontigname}{$depth}++;
		}
		if($blurb =~ /PROBLEM_INSERTION_OVERLAP/){
			%tempins = ();
			@temp_depth = ();
			#INS_671	13548_1_29	con1	Length:4619	AlignQuality:97.1	13548_1_29.1	326	GGGACTCGAAAGCGAAAGTT	con11	Length:1149	AlignQuality:25.2	13548_1_29.11	219	GTT	con3	Length:3014	AlignQuality:9.4	13548_1_29.3	481	GGGACTCGAAAGCGAAAGTT	PROBLEM_INSERTION_OVERLAP
			while($blurb =~ /con\S+\tLength:(\d+)\tAlignQuality:(\S+)\t(\S+)\t(\d+)\t(\S+)\t/g){
# UNUSED			$contiglength = $1;
				$alignquality = $2;
				$contig = $3;
				$pos = $4;
				$insert_seq = $5;
				$depth = $depth_hash{$contig}{$pos} // 0;
				# rank depth better than alignquality here - mainly because by definition longer insertions may be more troublesome to align and give lower AQ
				# even though the depth might be a lot higher ... 
				$tempins{$contig}{$depth}{$alignquality}{$insert_seq}++;
				if(exists($no_insertion_depth{$contig})){
					delete $no_insertion_depth{$contig};
				}
				push(@temp_depth,$depth);
			}
			
			if(scalar keys %no_insertion_depth){
				foreach $contig (keys %no_insertion_depth){
					foreach $depth (keys %{$no_insertion_depth{$contig}}){
						if($depth > $max_no_insertion_depth){
						$max_no_insertion_depth = $depth;
						}
					}
				}
			}

			@temp_depth_sorted = reverse sort {$a <=> $b} (@temp_depth);
			# i.e. if the depth of at least one sequence with an insertion is bigger than all those without carry on:
			if($max_no_insertion_depth < $temp_depth_sorted[0]){
				foreach $contig (keys %tempins){
				@tempdepths = reverse sort {$a <=> $b} (keys %{$tempins{$contig}});
				$maxdepth = $tempdepths[0];
				@tempaq =  reverse sort {$a <=> $b} (keys %{$tempins{$contig}{$maxdepth}});
				$maxaq = $tempaq[0];
					foreach $insert_seq (keys %{$tempins{$contig}{$maxdepth}{$maxaq}}){
					$probleminsertion{$insert_seq} += $depth;
					}
				}
			$fullmaxdepth = 0;
				foreach $insert_seq (keys %probleminsertion){
					if($probleminsertion{$insert_seq} > $fullmaxdepth){
					$fullmaxdepth = $probleminsertion{$insert_seq};
					}
				}
			@probleminsertion = ();
				foreach $insert_seq (keys %probleminsertion){
					if($probleminsertion{$insert_seq} == $fullmaxdepth){
					push(@probleminsertion,$insert_seq);
					}
				}
				if(scalar @probleminsertion == 1){
				$probleminsertion = $probleminsertion[0];
				$problemdepth = $fullmaxdepth;
				}
				else{
				$i = int rand(scalar @probleminsertion);
				$probleminsertion = $probleminsertion[$i];
				$problemdepth = $fullmaxdepth;
				}
			$sampleseq{$sample} .= $probleminsertion;
			
				for($i=0;$i<length($probleminsertion);$i++){
				$samplepos++;
				$samplemap{$sample}{$samplepos} = "INS_AFTER_".$referencepos;
				}
			
			} # end if max_no_insertion_depth < maxdepth of insertion...
			
			
			
			
		} #  end if problemoverlap and insertion
# so now just left with insertions...		
#INS_6776        12559_1_12      con1    Length:9016     AlignQuality:99.9       12559_1_12.1    6283    AATGATAGT
#INS_6237	10065_1_1	con2	Length:3459	AlignQuality:100	10065_1_1.2	272	GGATCAGGA
		elsif($blurb =~ /con\S+\tLength:(\d+)\tAlignQuality:(\S+)\t(\S+)\t(\d+)\t(\S+)/){
		$alignquality = $2;
		$contig = $3;
		$pos = $4;
		$insert_seq = $5;
		$insert_depth = $depth_hash{$contig}{$pos} // 0;
		$tempins{$contig}{$depth}{$alignquality}{$insert_seq}++;
				
			if(exists($no_insertion_depth{$contig})){
			delete $no_insertion_depth{$contig};
			}

			if(scalar keys %no_insertion_depth){
				foreach $contig (keys %no_insertion_depth){
					foreach $depth (keys %{$no_insertion_depth{$contig}}){
						if($depth > $max_no_insertion_depth){
						$max_no_insertion_depth = $depth;
						}
					}
				}
			}
			if($max_no_insertion_depth < $insert_depth){
			
			$sampleseq{$sample} .= $insert_seq;
				for($i = 0; $i < length($insert_seq); $i++){
				$samplepos++;
				$samplemap{$sample}{$samplepos} = "INS_AFTER_".$referencepos;
				}
			}
		}
	} # end elsif INS_ 
}

close(LASTZ_ANALYSED) || die $!;

open(OUT,">$output") || die $!;
foreach $sample (sort keys %sampleseq){
print OUT ">".$sample.".genome\n",$sampleseq{$sample},"\n";
}
close(OUT) || die $!;

$samplemap = $output.".samplemap.txt";
open(MAP,">$samplemap") || die $!;
foreach $sample (sort keys %samplemap){
	foreach $samplepos (sort {$a <=> $b} keys %{$samplemap{$sample}}){
	print MAP "$sample\t$samplepos\t".$samplemap{$sample}{$samplepos}."\n";
	}
}
close(MAP) || die $!;
