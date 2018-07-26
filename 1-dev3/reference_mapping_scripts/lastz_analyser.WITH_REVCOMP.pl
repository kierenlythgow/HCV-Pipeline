#!/usr/bin/perl -s


## using switches now to make pipeline easier to work with.

## i.e. perl -w -s lastz_analyser.pl
## -reference_fasta_file=reference.fasta
## -sample_fasta_file=sample.fasta
## -lastz_results_file=sample.lastz
## -log_file=logs/sample.lastz_analyser.log
## -cutoff=50000
## -with_revcomp=yes or no
## -output=sample.lastz_analysed_file
## Use cutoff to remove low-scoring lastz alignments
## This prevents secondary alignments from repeat regions to front/back of reference genome.

##### Aim is to map reference to each of the contigs so can compare differences at each
##### reference position where there may be contig overlap and also to identify gaps in the alignment
##### that perhaps we can fill in using reference-based bwa/smalt later on, where IVA assembly may have failed.
##### e.g. ref 530 -> contig_1 433 -> A
#####      ref 530 -> contig_2 204 -> C

##### There may be cases of dual infection where we want to use two reference sequences.  This shouldn't be a problem - simply run the script twice using
##### a different reference each time.  Then we'll have maps to both references that will (hopefully) mean we can generate two genomes down the line to
##### remap to using smalt/bwa simultaneously.

if(!$with_revcomp){
	$with_revcomp="no";
}
open(OUT,">$output") || die $!;


# Can change this is we're using one log file per sample per pipeline run.
# Currently I'm keeping this logfile separate for ease of validation
open(LOG,">$log_file") || die $!;

# First load the reference sequence in.
open(REF,"$reference_fasta_file") || die $!;
while(<REF>){
	if(!/>/){
		chomp($_);
		$reference_sequence .= $_;
	}
}
close(REF) || die $!;

# Now load the contigs in.
open(CONTIG_SEQS,"$sample_fasta_file") || die $!;
while(<CONTIG_SEQS>){
	chomp($_);
	if(/>/){
		$contigname = $_;
		$contigname =~ s/\>//;

		#Following steps can be removed once the read files, contig files and contig names all follow a standard notation.  WTSI working on solving this (Jan 2015).
		#Next step to sort out Swee's rather unusual contig-naming convention for something like 10% of the BEEHIVE samples
		$contigname =~ s/EGAR\d+\.//;
		#And this to remove any remaining annotation from the IVA0.7 vs. IVA0.5 comparison
		$contigname =~ s/IVA07__//;
	}
	else{
		$contig_sequences{$contigname}.=$_;
	}
}
close(CONTIG_SEQS) || die $!;


# Now take in full lastz output i.e. the form that gives the individual alignments per contig per sample
$alignment_part_number = 0;

#while(<>){
open(LASTZ,"$lastz_results_file") || die $!;
while(<LASTZ>){

#this will need modding depending on where and how lastz was run
#lastz.v1.03.54 hxb2.newtopline.fasta PG14-UG002011-S00011.fa
#lastz.v1.03.54 ../hxb2.newtopline.fasta ../10065_1_1.fa
### FOR THE PIPELINE WE SHOULD HAVE ONE LASTZ OUTPUT PER SAMPLE, RATHER THAN A MULTIPLE SAMPLES WHICH WE GENERATED FOR BATCH JOBS BEFORE

### i.e. in the form lastz.v1.03.54 $reference.fasta $sample.fasta
	#if(/\"lastz\S+.*\s+.*\/+(\S+)\.fa/){
	if(/\"lastz\S+.*\s+(\S+)\.fa/){
		$sample = $1;
		# Remove the full path before $sample
		# i.e. /home/data/12345_1_10.fa becomes 12345_1_10
		$sample =~ s/.*\///;
	#this added when contig files were combined and the filename ended up as "sample.contigs.fasta"
		$sample =~ s/\.contigs//;
		# Remove any trailing numbers
		# i.e. /home/data/12345_1_10.1.fa becomes 12345_1_10 after these two steps
		$sample =~ s/\.\d+$//;
		$alignment_number = 0;
		$alignment_on = 0;
		$samplenamehash{$sample}++;
	}
	elsif(/\s+\"\>(\S+)\.(\d+)/){
# # # # # # Slight problem in analysis here - we noticed IVA very occasionally fails in assembly and creates reverse complement alignments
# # # # # # AAAAAAAAGHHHHHH ">.13591_1_37.5 (reverse complement)" FAILS SO USES PRIOR CONTIG NAME UGH
# # # # # # Can adapt script later if required to make use of reverse complement alignments - if anything they're a sign that assembly has failed
# # # # # # so perhaps leave out and the lack of sequence later should flag the sample up as being a bit iffy...
		$contignumber = $2;
		if($_ =~/reverse\ complement/){
	    		$rc = "_RC";
			$revcomp = 1;
			if($with_revcomp eq "yes"){
				$revcomp=0;
			}
		}
		else{
			$rc = "";
			$revcomp = 0;
		}
		$contigname = $sample.'.'.$contignumber;
		$contig_sequence = $contig_sequences{$contigname};
		if($rc ne ""){
			$contig_sequence = reverse $contig_sequence;
			$contig_sequence =~ tr/ACGT/TGCA/;
			$contigname .= "_RC";
			$contig_sequences{$contigname} = $contig_sequence;
		}
	}
	### next line is a score filter to remove the lastz alignments with low scores that tend to be secondary alignments to repeat regions...
	elsif(/\s+s\s+(\d+)/){
		$score = $1;
	}
	elsif(/\s+b\s+(\d+)\s+(\d+)/){
		$reference_full_alignment_beginning = $1;
		$contig_full_alignment_beginning = $2;
		$alignment_number++;
		$contig_part_alignment_seq = "";
		$reference_part_alignment_seq = "";
		$alignment_on = 1;
	}
	elsif(/\s+e\s+(\d+)\s+(\d+)/){
		$reference_full_alignment_end = $1;
		$contig_full_alignment_end = $2;
		$alignment_part_number = 0;
# recalc'd alignment quality to use contig beginning and end instead of reference as makes more sense and can't now get qualities>100%!
		$alignment_quality{$sample}{$contigname}{$alignment_number} = 100*sprintf("%.3f", ($contig_full_alignment_end - $contig_full_alignment_beginning + 1) / length($contig_sequences{$contigname}));
		if(($score > $cutoff) && $revcomp == 0){
			for($referencepos = $reference_full_alignment_beginning; $referencepos<= $reference_full_alignment_end; $referencepos++){
				$hyphen_hash{$sample}{$referencepos} = '-';
			}
		}
	}
	elsif(/\s+l\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)/  && ($score > $cutoff) && $revcomp == 0){
		if($alignment_part_number>0){
		$old_contig_part_stop = $contig_part_stop;
		$old_reference_part_stop = $reference_part_stop;
		}
		$alignment_part_number++;
		$reference_part_start = $1;
		$contig_part_start = $2;
		$reference_part_stop = $3;
		$contig_part_stop = $4;
		$alignment_on = 1;
		for($i=$reference_part_start;$i<=$reference_part_stop;$i++){
			$temp_contig_pos = $contig_part_start+$i-$reference_part_start;
			#store map from reference pos to contig_pos
			$mapping{$i}{$sample}{$contigname}{$alignment_number} = $temp_contig_pos;
			$temp_contig_pos++;
		}
		if($alignment_part_number>1){
			if(($contig_part_start-$old_contig_part_stop)>1){
				for($i=1;$i<($contig_part_start-$old_contig_part_stop);$i++){
					$reference_part_alignment_seq .= '.';
				}
				$contig_part_alignment_seq .= substr($contig_sequence,$old_contig_part_stop,($contig_part_start-$old_contig_part_stop-1));
				$insertionseq{$sample}{$contigname}{$alignment_number}{$old_reference_part_stop} = substr($contig_sequence,$old_contig_part_stop,($contig_part_start-$old_contig_part_stop-1));
				$insertionposLHS{$sample}{$contigname}{$alignment_number}{$old_reference_part_stop} = $old_contig_part_stop+1;
			}
			if(($reference_part_start-$old_reference_part_stop)>1){
# UNUSED			$temp_del_seq = "";
				for($i=1;$i<($reference_part_start-$old_reference_part_stop);$i++){
					$contig_part_alignment_seq .= '.';
					$reference_part_alignment_seq .= substr($reference_sequence,($old_reference_part_stop+$i-1),1);
					$deletionseq{$sample}{$contigname}{$alignment_number} {$old_reference_part_stop+$i}  = substr($reference_sequence,($old_reference_part_stop+$i-1),1);
					$deletionposLHS{$sample}{$contigname}{$alignment_number} {$old_reference_part_stop+$i} = $old_contig_part_stop;
				}
				$deletionseqfull{$sample}{$contigname}{$alignment_number}{$old_reference_part_stop+1} = substr($reference_sequence,$old_reference_part_stop,($reference_part_start-$old_reference_part_stop-1));
				$del_start_to_end{$sample}{$contigname}{$alignment_number}{$old_reference_part_stop+1} = ($reference_part_start - 1);
			}
		}
		$contig_part_alignment_seq .= substr($contig_sequence,($contig_part_start-1),($contig_part_stop-$contig_part_start+1));
		$reference_part_alignment_seq .= substr($reference_sequence,($reference_part_start-1),($reference_part_stop-$reference_part_start+1));
	}
	elsif(/\}/ && ($alignment_on==1)){
		#print OUT output then tell script alignment parts are OFF
		print LOG "$sample\t$contigname\t$alignment_number\nreferencealign\t$reference_part_alignment_seq\ncontigalign\t$contig_part_alignment_seq\n";
		$alignment_on = 0;
	}
}

close(LASTZ) || die $!;

foreach $sample (keys %insertionposLHS){
	foreach $contig (sort keys %{$insertionposLHS{$sample}}){
		foreach $alignment (sort {$a <=> $b} keys %{$insertionposLHS{$sample}{$contig}}){
			foreach $pos (sort {$a <=> $b} keys %{$insertionposLHS{$sample}{$contig}{$alignment}}){
				print LOG "INSERTION:$sample;$contig;referenceposLHS:$pos;contigposLHS:".$insertionposLHS{$sample}{$contig}{$alignment}{$pos}.";sequence:".$insertionseq{$sample}{$contig}{$alignment}{$pos}."\n";
			}
		}
	}
}

foreach $sample (keys %deletionposLHS){
	foreach $contig (sort keys %{$deletionposLHS{$sample}}){
		foreach $alignment (sort {$a <=> $b} keys %{$deletionposLHS{$sample}{$contig}}){
			foreach $pos (sort {$a <=> $b} keys %{$deletionposLHS{$sample}{$contig}{$alignment}}){
				print LOG "DELETION:$sample;$contig;referencepos:$pos;contigposLHS:".$deletionposLHS{$sample}{$contig}{$alignment}{$pos}.";sequence:".$deletionseq{$sample}{$contig}{$alignment}{$pos}."\n";
			}
			foreach $startpos (sort {$a <=> $b} keys %{$deletionseqfull{$sample}{$contig}{$alignment}}){
				$endpos = $del_start_to_end{$sample}{$contig}{$alignment}{$startpos};
				print LOG "FULL_DELETION:$sample;$contig;referencestartpos:$startpos;referenceendpos:$endpos;sequence_full:".$deletionseqfull{$sample}{$contig}{$alignment}{$startpos}."\n";
			}
		}
	}
}
###Old line: when we were just using HXB2 as the reference, which (perhaps unsurprisingly given the code) is 9719nt long...
###@referencepos = (1..9719);
@referencepos = (1..length($reference_sequence));
foreach $sample(sort keys %samplenamehash){
	foreach $referencepos (@referencepos){
		$referencenucl = substr($reference_sequence,($referencepos-1),1);
		print OUT "$referencepos\t".$referencenucl;
		print OUT "\t$sample";
		%tempnucl = ();
		$missingdataflag = 0;
		foreach $contig (sort keys %{$mapping{$referencepos}{$sample}}){
			if($contig =~ /.*\.(\d+)$/){
				$contignumber = "con$1";
			}
			elsif($contig =~ /.*\.(\d+)_RC$/){
				$contignumber = "con$1"."_RC";
			}
			$missingdataflag = 1;
			$contiglength = length($contig_sequences{$contig});
			foreach $alignment (sort {$a<=>$b} keys %{$mapping{$referencepos}{$sample}{$contig}}){
				$contigpos = $mapping{$referencepos}{$sample}{$contig}{$alignment};
				$contignucl = substr($contig_sequences{$contig},($contigpos-1),1);
				$tempnucl{$contignucl}++;
				print OUT "\t$contignumber\tLength:$contiglength\tAlignQuality:".$alignment_quality{$sample}{$contig}{$alignment};
				print OUT "\t$contig\t$contigpos\t$contignucl";
			}
		}
		if(scalar(keys %tempnucl)>1){
			print OUT "\tPROBLEM_OVERLAP\n";
		}
		elsif(scalar(keys %tempnucl)==1){
			@tempres = keys %tempnucl;
			if($tempres[0] ne $referencenucl){
				print OUT "\tMUTATION:$referencenucl".$referencepos.$tempres[0]."\n";
			}
			else{
				print OUT "\n";
			}
		}
		elsif($missingdataflag == 0){
			if(exists($hyphen_hash{$sample}{$referencepos})){
				print OUT "\tDELETION_IN_ALIGNMENT_TO_REFERENCE\n";
			}
			else{
				print OUT "\tNO_ALIGNMENT_TO_REFERENCE_FOR_THIS_POSITION\n";
			}
		}
		else{
			print OUT "\n";
		}

		#NOW TO DEAL WITH INSERTIONS:
		$line = "INS_$referencepos\t$sample";
		$restofline = "";
		%tempins = ();
		foreach $contig (sort keys %{$mapping{$referencepos}{$sample}}){
			if($contig =~ /.*\.(\d+)$/){
				$contignumber = "con$1";
			}
			$contiglength = length($contig_sequences{$contig});
			foreach $alignment (sort {$a<=>$b} keys %{$mapping{$referencepos}{$sample}{$contig}}){
				$contigpos = $mapping{$referencepos}{$sample}{$contig}{$alignment};
				$contignucl = substr($contig_sequences{$contig},($contigpos-1),1);
				$tempnucl{$contignucl}++;
				if(exists($insertionposLHS{$sample}{$contig}{$alignment}{$referencepos})){
					$insertionseq = $insertionseq{$sample}{$contig}{$alignment}{$referencepos};
					$restofline .= "\t$contignumber\tLength:$contiglength\tAlignQuality:".$alignment_quality{$sample}{$contig}{$alignment};
					$restofline .= "\t$contig\t$contigpos\t";
					$restofline .= $insertionseq;
					$tempins{$insertionseq}++;
				}
			}

		}
		if($restofline ne ""){
			print OUT $line,$restofline;
			if(scalar(keys %tempins)>1){
				print OUT "\tPROBLEM_INSERTION_OVERLAP\n";
			}
			else{
				print OUT "\n";
			}
		}
	}
}

close(LOG) || die $!;
close(OUT) || die $!;

