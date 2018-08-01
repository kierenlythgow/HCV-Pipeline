#!/usr/bin/perl -s

#usage for first round of consensus/mv generation:
# perl -w cons_mv.pl 
# -mpileup=13774_1_54.genome.fasta.smalt.sorted.mpileup 
# -reference_fasta=13774_1_54.genome.fasta 
# -consensus_out=13774_1_54.consensus1.preNcut.fasta 
# -mv_out=13774_1_54.genome.fasta.mv 
# -base_freq_out=13774_1_54.genome.fasta.basefreqs.tsv

# later on use:
# perl -w cons_mv.pl 
# -mpileup=13774_1_54.consensus1.fasta.smalt.sorted.mpileup 
# -reference_fasta=13774_1_54.consensus1.fasta 
# -consensus_out=13774_1_54.consensus2.preNcut.fasta 
# -mv_out=13774_1_54.consensus1.fasta.mv 
# -base_freq_out=13774_1_54.consensus1.fasta.basefreqs.tsv

# i.e. for each iteration the final consensus out file X goes up by one...
# if there are no majority errors in sample.consensus(X).fasta.mv, then the sample.consensus(X+1).fasta and sample.consensus(X).fasta 
# sequences *should* be identical! 
#
# That said...
# If a base is called at a depth below the threshold, it will get reported as an N.
# Once a run of N's get excised, occasionally reads map differently or "across the divide" 
# especially at low, contaminant depth levels, so an additional round of mapping may
# throw up a consensus error after a round where there are none...
# Solution: throw in a flag to report A/C/G/T -> N as a majority error
# Problem: low depth sequences can end up shrinking until you end up with a run of N's
# Potential counter: use a cut-off of 1 until the very last step.  *Should* work.

if(!$mpileup){
die "Error: script needs an mpileup file\nUsage: perl -w -s cons_mv.pl -mpileup=PATH_TO_MPILEUP_FILE/MPILEUP_FILENAME\n\n";
}
if(!$reference_fasta){
die "Error: script needs a reference FASTA file\nUsage: perl -w -s cons_mv.pl -reference_fasta=PATH_TO_REFERENCE_FASTA_FILE/REFERENCE_FASTA_FILENAME\n\n";
}
if(!$mv_freq_cutoff){
$mv_freq_cutoff = 0.05;
}
if(!$mv_overall_depth_cutoff){
$mv_overall_depth_cutoff = 1000;
}
if(!$mv_variant_depth_cutoff){
$mv_variant_depth_cutoff = 50;
}
if(!$cons_depth_cutoff){
$cons_depth_cutoff = 100;
}
if(!$sliding_window_size){
$sliding_window_size = 300;
}


## $consensus_out : FASTA output of consensus sequence
## $mv_out: flatfile record of minority variants and consensus errors

if(!$consensus_out){
$consensus_out = $mpileup;
$consensus_out =~ s/mpileup/consensus\.fasta/;
}
if(!$mv_out){
$mv_out = $mpileup;
$mv_out =~ s/mpileup/mv/;
}
if(!$base_freq_out){
$base_freq_out = $mv_out;
$base_freq_out =~ s/\.mv$/\.basefreqs.tsv/;
}



$logcons = $consensus_out.".logfile";

$|++;

open(REF,$reference_fasta) || die $!;
while(<REF>){
	if(/>(.*)\n/){
	$refname = $1;
	}
	else{
	chomp($_);
	$reference_sequence .= $_;
	}
}
close(REF) || die $!;

# store reference bases and initialise separate hash with depth set to 0 - this will get updated with depth of coverage when parsing mpileup file later.
# can then loop through this at the end to work out which bases need correcting and which are not covered by any or sufficient numbers of reads to warrant inclusion in final genome...
for($i = 0; $i < length($reference_sequence); $i++){
$consensus_sequence{($i+1)} = substr($reference_sequence,$i,1);
$consensus_depth{($i+1)}  = 0;
}

open(MV,">$mv_out") || die $!;
open(CONS,">$consensus_out") || die $!;

open(IN,"$mpileup")|| die $!;
print MV "Sample\tSampleConsensusPosition\tConsensusBase\tDepth\tType\tVariantBase\tVariantDepth\tVariantFrequency\tVariantFwdRevRatio\tConsensusFwdRevRatio\n";

%consensus_ins = ();
%consensus_del = ();
%acgtd = ();

while(<IN>){
#94C_HIV_T146_Cons       400     A       2044    ,-1g,,,,,,,,,,-1g,,-1g,-1g,,,,-1g,,,,,,,,,,-1g,,

@line = split(/\t/,$_);
$seq = $line[0];
$pos = $line[1];
$refbase = $line[2];
$depth = $line[3];
$consensus_depth{$pos} = $depth;
$baselist = $line[4];
%tempindel = ();
%tempvar = ();
%minindel = ();
%minindelseq = ();
$identical = 0;
$multidel = 0;
%fwd = ();
%rev = ();
## Initialise the fwd and rev count hash for bases so that if for some reason there are no reads that match the refbase, then you don't end up with 
## uninitialised fwd or rev counts when calculating the fwd/(fwd+rev) ratio later on...
	
	foreach $base (qw(A C G T B D H I K M N R S U V W X Y)){
	$fwd{$base} = 0;
	$rev{$base} = 0;
	$fwd{lc($base)} = 0;
	$rev{lc($base)} = 0;
	}
	
	
	
#don't care about strand - this upper case command makes regular expression twice as efficient
#and a lot easier to read
$string = uc($baselist);
$string =~ s/\^\S//g;
$string =~ s/\$//g;
$identical = 0;
$newbaselist = $baselist;
$newbaselist =~ s/\^\S//g;
$newbaselist =~ s/\$//g;
	while($string =~ /(,|\.)/g){
	$identical++;

### Argh, made all stored refbases here upper case - only an issue now we're using (from May2016) lower
### case bases in sequences to denote iffy read depth...
		if($1 =~ /\./){
		$fwd{uc($refbase)}++;
		}
		else{
		$rev{uc($refbase)}++;
		}
	$acgtd{$pos}{uc($refbase)}++;
	}

	while($string =~ /(\+|\-)(\d+(A|C|G|T|N)+)/g){

##  need to edit +1TTTG -2CACG etc. as the TTG and CG are substitutions after insertions/deletions...
	$plusminus = $1;
	$indelsub = $2;
		if($indelsub =~ /^(\d+)(\S+)/){
		$indellength = $1;
		$seq_indelsubs = $2;
		$actualindel = substr($seq_indelsubs,0,$indellength);
		}
		else{
		print MV "BBBBBBBBBBBBB___MALFORMED_CIGAR_CODE:";
		}
	$tempindel{$plusminus}{"$indellength"."$actualindel"}++;
		for($i = 1; $i <= $indellength ; $i++){
		$minindel{$plusminus}{$i}++;
		$tempbase = substr($actualindel,($i-1),1);
		$minindelseq{$plusminus}{$i}{$tempbase}++;
		}
	$remainder = substr($indelsub,($indellength+length($indellength)));
		while($remainder =~ /(A|C|G|T|N)/g){
		$tempvar{$1}++;
		}
	}

	
if(($depth >= $mv_overall_depth_cutoff) || ($depth > 2*$identical)){

#remove indels so should just be left with substitutions AND multidel "*"'s
$string =~ s/(\+|\-)(\d+)((A|C|G|T|N)+)//g;
	while($string =~ /(A|C|G|T|N)/g){
	$tempvar{$1}++;
	}

	while($string =~ /\*/g){
	$multidel++;
	}

###NOW LOOP THROUGH tempvar and find the ACG or T with the highest frequency. if this > $identical, reset refbase to this max and call a maj error below...
$tempmax = $identical;
$temprefbase = $refbase;
	foreach $var (sort keys %tempvar){
		if($tempvar{$var} > $tempmax){
		$tempmax = $tempvar{$var};
		$temprefbase = $var;
		}
	}
$majerror = "";
###added to cases where there's a majority deletion - otherwise wrongrefbase never gets set. tidy below later.
$wrongrefbase=$refbase;
	if($temprefbase ne $refbase){
	$wrongrefbase = $refbase;
	$refbase = $temprefbase;
	$majerror = "MajError";
	}
$consensus_sequence{$pos} = $refbase;


####extra loops at the start to work out fwd and reverse counts
####really should rewrite code to not uppercase the string until when it's actually needed - but time is of the essence - tidy that up later!!!



	foreach $var (sort keys %{$tempindel{'+'}}){
	$fwd{$var} = 0;
	$fwd_var = '\+'.$var;	
		while($newbaselist =~ /$fwd_var/g){
		$fwd{$var}++;
		$newbaselist = $`.$';
		}
	$rev{$var} = 0;	
	$rev_var = '\+'.lc($var);
		while($newbaselist =~ /$rev_var/g){
		$rev{$var}++;
		$newbaselist = $`.$';
		}
	}
	
	foreach $var (sort keys %{$tempindel{'-'}}){
	$fwd_var = '\-'.$var;	
	$fwd{$var} = 0;
		while($newbaselist =~ /$fwd_var/g){
		$fwd{$var}++;
		$newbaselist = $`.$';
		}
	$rev{$var} = 0;	
	$rev_var = '\-'.lc($var);
		while($newbaselist =~ /$rev_var/g){
		$rev{$var}++;
		$newbaselist = $`.$';
		}
	}
	foreach $var (sort keys %tempvar){
	$fwd_var = uc($var);
	$fwd{$var} = 0;
		while($newbaselist =~ /$fwd_var/g){
		$fwd{$var}++;
		}
	$rev{$var} = 0;	
	$rev_var = lc($var);
		while($newbaselist =~ /$rev_var/g){
		$rev{$var}++;
		}
	}
	
$total_fwd_refbase = $fwd{$refbase};	
$total_rev_refbase = $rev{$refbase};	
	if(($total_rev_refbase+$total_fwd_refbase) > 0){
	$fr_ratio_refbase = sprintf("%.2f",($total_fwd_refbase/($total_rev_refbase+$total_fwd_refbase)));
	}
	else{
	$fr_ratio_refbase = 0;
	}

###Now go through substitution variants, stored in tempvar...	
	
	foreach $var (sort keys %tempvar){

	$freq = $tempvar{$var} / $depth; 
	$rounded = sprintf("%.3f", $freq);
		if($rounded >= $mv_freq_cutoff){
			if($tempvar{$var} >= $mv_variant_depth_cutoff){
			$fr_ratio_mv = sprintf("%.2f",($fwd{$var}/($fwd{$var}+$rev{$var})));
			print MV "$seq\t$pos\t$refbase\t$depth\t".$majerror."Sub\t$var\t$tempvar{$var}\t$rounded\t$fr_ratio_mv\t$fr_ratio_refbase\n";
			}
			if($freq > $identical/$depth){

		#	print MV out new majority base in refbase position and refbase position in minority variant position plus the number of times reference base seen (i.e. .'s and ,'s summed)
		#	need this line to store new min var stats as %incorrect ref != 1-%majority variant necessarily does it...
				### need this if in cases where overall depth and variant depth are low but consensus is wrong.
				### e.g. depth 95, majority variant=48, min var=47
				### as the variant won't have printed out above...
				if($tempvar{$var} < $mv_variant_depth_cutoff){
				$fr_ratio_mv = sprintf("%.2f",($fwd{$var}/($fwd{$var}+$rev{$var})));
				print MV "$seq\t$pos\t$refbase\t$depth\t".$majerror."SubLow\t$var\t$tempvar{$var}\t$rounded\t$fr_ratio_mv\t$fr_ratio_refbase\n";
				$rounded2 = sprintf("%.3f", ($identical/$depth));
### CORRECTING ERROR WHERE WRONGREFBASE HAS NO COVERAGE IN EITHER DIRECTION
					if(($rev{$wrongrefbase}+$fwd{$wrongrefbase}) != 0){
					$fr_ratio_wrongrefbase = sprintf("%.2f",($fwd{$wrongrefbase}/($rev{$wrongrefbase}+$fwd{$wrongrefbase})));
					}
					else{
					$fr_ratio_wrongrefbase = 0;
					}
				print MV "$seq\t$pos\t$refbase\t$depth\tMajErrorSubLow\t$wrongrefbase\t$identical\t$rounded2\t$fr_ratio_wrongrefbase\t$fr_ratio_refbase\n";
				}
				else{
				
				$rounded2 = sprintf("%.3f", ($identical/$depth));
### CORRECTING ERROR WHERE WRONGREFBASE HAS NO COVERAGE IN EITHER DIRECTION
					if(($rev{$wrongrefbase}+$fwd{$wrongrefbase}) != 0){
					$fr_ratio_wrongrefbase = sprintf("%.2f",($fwd{$wrongrefbase}/($rev{$wrongrefbase}+$fwd{$wrongrefbase})));
					}
					else{
					$fr_ratio_wrongrefbase = 0;
					}
				print MV "$seq\t$pos\t$refbase\t$depth\tMajErrorSub\t$wrongrefbase\t$identical\t$rounded2\t$fr_ratio_wrongrefbase\t$fr_ratio_refbase\n";
				}
			}
		}
	$acgtd{$pos}{$var} = $tempvar{$var};
	}	
		
	
	foreach $var (sort keys %{$tempindel{'+'}}){

        $freq = $tempindel{'+'}{$var} / $depth; 
        $rounded = sprintf("%.3f", $freq);
	        if($rounded >= $mv_freq_cutoff){
			if($tempindel{'+'}{$var} >= $mv_variant_depth_cutoff){
			$fr_ratio_mv = sprintf("%.2f",($fwd{$var}/($fwd{$var}+$rev{$var})));
			print MV "$seq\t$pos\t$refbase\t$depth\t".$majerror."Ins\t$var\t$tempindel{'+'}{$var}\t$rounded\t$fr_ratio_mv\t$fr_ratio_refbase\n";
			}
			if($freq >= 0.5){
		#	print MV out new majority base in refbase position and refbase position in minority variant position plus the number of times reference base seen (i.e. .'s and ,'s summed)
		#	need this line to store new min var stats as %incorrect ref != 1-%majority variant necessarily does it...
				### need this if in cases where overall depth and variant depth are low but consensus is wrong.
				### e.g. depth 95, majority variant=48, min var=47
				### as the variant won't have printed out above...
				if($tempindel{'+'}{$var} < $mv_variant_depth_cutoff){
				$fr_ratio_mv = sprintf("%.2f",($fwd{$var}/($fwd{$var}+$rev{$var})));
				print MV "$seq\t$pos\t$refbase\t$depth\tInsLow\t$var\t$tempindel{'+'}{$var}\t$rounded\t$fr_ratio_mv\t$fr_ratio_refbase\n";
				$rounded2 = sprintf("%.3f", ($identical/$depth));
### CORRECTING ERROR WHERE WRONGREFBASE HAS NO COVERAGE IN EITHER DIRECTION
					if(($rev{$wrongrefbase}+$fwd{$wrongrefbase}) != 0){
					$fr_ratio_wrongrefbase = sprintf("%.2f",($fwd{$wrongrefbase}/($rev{$wrongrefbase}+$fwd{$wrongrefbase})));
					}
					else{
					$fr_ratio_wrongrefbase = 0;
					}
				print MV "$seq\t$pos\t$var\t$depth\tMajErrorInsLow\t$refbase\t$identical\t$rounded2\t$fr_ratio_wrongrefbase\t$fr_ratio_refbase\n";
				}
				else{
				$rounded2 = sprintf("%.3f", ($identical/$depth));
			#	print LOG "ins id $identical, $rounded, $rounded2\n";
### CORRECTING ERROR WHERE WRONGREFBASE HAS NO COVERAGE IN EITHER DIRECTION
					if(($rev{$wrongrefbase}+$fwd{$wrongrefbase}) != 0){
					$fr_ratio_wrongrefbase = sprintf("%.2f",($fwd{$wrongrefbase}/($rev{$wrongrefbase}+$fwd{$wrongrefbase})));
					}
					else{
					$fr_ratio_wrongrefbase = 0;
					}
				print MV "$seq\t$pos\t$var\t$depth\tMajErrorIns\t$refbase\t$identical\t$rounded2\t$fr_ratio_wrongrefbase\t$fr_ratio_refbase\n";
				}
			$consensus_ins{$pos} = $var;
			}
		}
	}
	
	
	foreach $var (sort keys %{$tempindel{'-'}}){
	
        $freq = $tempindel{'-'}{$var} / $depth;
        $rounded = sprintf("%.3f", $freq);
		if($rounded >= $mv_freq_cutoff){
			if($tempindel{'-'}{$var} >= $mv_variant_depth_cutoff){
			$fr_ratio_mv = sprintf("%.2f",($fwd{$var}/($fwd{$var}+$rev{$var})));
			$rounded2 = sprintf("%.3f", ($identical/$depth));
			print MV "$seq\t$pos\t$refbase\t$depth\t".$majerror."Del\t$var\t$tempindel{'-'}{$var}\t$rounded\t$fr_ratio_mv\t$fr_ratio_refbase\n";
			}
			if($freq >= 0.5){
		#	print MV out new majority base in refbase position and refbase position in minority variant position plus the number of times reference base seen (i.e. .'s and ,'s summed)
		#	need this line to store new min var stats as %incorrect ref != 1-%majority variant necessarily does it...
			### need this if in cases where overall depth and variant depth are low but consensus is wrong.
			### e.g. depth 95, majority variant=48, min var=47
			### as the variant won't have printed out above...
				if($tempindel{'-'}{$var} < $mv_variant_depth_cutoff){
				$fr_ratio_mv = sprintf("%.2f",($fwd{$var}/($fwd{$var}+$rev{$var})));
				print MV "$seq\t$pos\t$refbase\t$depth\tDelLow\t$var\t$tempindel{'-'}{$var}\t$rounded\t$fr_ratio_mv\t$fr_ratio_refbase\n";
				$rounded2 = sprintf("%.3f", ($identical/$depth));
### CORRECTING ERROR WHERE WRONGREFBASE HAS NO COVERAGE IN EITHER DIRECTION
					if(($rev{$wrongrefbase}+$fwd{$wrongrefbase}) != 0){
					$fr_ratio_wrongrefbase = sprintf("%.2f",($fwd{$wrongrefbase}/($rev{$wrongrefbase}+$fwd{$wrongrefbase})));
					}
					else{
					$fr_ratio_wrongrefbase = 0;
					}
				print MV "$seq\t$pos\t$var\t$depth\tMajErrorDelLow\t$refbase\t$identical\t$rounded2\t$fr_ratio_wrongrefbase\t$fr_ratio_refbase\n";
				}
				else{
				$rounded2 = sprintf("%.3f", ($identical/$depth));
### CORRECTING ERROR WHERE WRONGREFBASE HAS NO COVERAGE IN EITHER DIRECTION
					if(($rev{$wrongrefbase}+$fwd{$wrongrefbase}) != 0){
					$fr_ratio_wrongrefbase = sprintf("%.2f",($fwd{$wrongrefbase}/($rev{$wrongrefbase}+$fwd{$wrongrefbase})));
					}
					else{
					$fr_ratio_wrongrefbase = 0;
					}
				print MV "$seq\t$pos\t$var\t$depth\tMajErrorDel\t$refbase\t$identical\t$rounded2\t$fr_ratio_wrongrefbase\t$fr_ratio_refbase\n";
				}
			$consensus_del{$pos} = $var;
			}
		}
	}
### BELOW:  extra lines put in to account for the fact that you might have identity 45 bases, -1A:10, -2AT:40, -3ATC:25 in which case total depth = 120, 
### All deletions currently < identity so come up as MV's BUT freq of AT LEAST 1bp deletion = 10+40+25 = 75/120 = MajError... AT LEAST 2bp = 65/120=MajError also...
### But 3bp = 25/120 = MV only.  By storing the deletions earlier on in minindel, we can sort by increasing length of deletion and report the longest minimal deletion > identity.

### Keeping previous code as it might be interesting to know that of a 60% majerror deletion, 30% was -3ATC and 30% was -12ATCGGAGATGAC say - otherwise this info gets lost.
### NB in terms of adjusting the consensus, the final longest minimal deletion overwrites whatever was stored in consensus_del above.  Eventually will remove that line above
### but for now leaving it in so can switch code below on and off without changing otherwise fine code above.
### Then for each variant output an overall fwd/rev ratio, consensus ratio and variant ratio...
	$min_indel_seq = "";
	foreach $del_length (sort {$a <=> $b} keys %{$minindel{"-"}}){
	$minindelfreq = $minindel{"-"}{$del_length} / $depth;
	$minindelfreqrounded = sprintf("%.3f", $minindelfreq);
	#nb there should only be one key here - the reference base!
	$min_indel_seq .= (keys %{$minindelseq{"-"}{$del_length}})[0];
	$var = $del_length.$min_indel_seq;
		if($minindelfreqrounded >= $mv_freq_cutoff){
			if($minindel{"-"}{$del_length} >= $mv_variant_depth_cutoff){
			$identityrounded2 = sprintf("%.3f", ($identical/$depth));
			print MV "$seq\t$pos\t$refbase\t$depth\t".$majerror."MinDel\t$var\t".$minindel{"-"}{$del_length}."\t$minindelfreqrounded\t0\t$fr_ratio_refbase\n";
			}
			if($minindelfreq >= 0.5){
				if($minindel{"-"}{$del_length} < $mv_variant_depth_cutoff){
				print MV "$seq\t$pos\t$refbase\t$depth\tMinDelLow\t$var\t".$minindel{"-"}{$del_length}."\t$minindelfreqrounded\t0\t$fr_ratio_refbase\n";
				$identityrounded2 = sprintf("%.3f", ($identical/$depth));
				print MV "$seq\t$pos\t$var\t$depth\tMajErrorMinDelLow\t$refbase\t$identical\t$identityrounded2\t0\t$fr_ratio_refbase\n";
				}
				else{
				$identityrounded2 = sprintf("%.3f", ($identical/$depth));
				print MV "$seq\t$pos\t$var\t$depth\tMajErrorMinDel\t$refbase\t$identical\t$identityrounded2\t0\t$fr_ratio_refbase\n";
				}
			$consensus_del{$pos} = $var;
			}
		}
	}
### NOW DO SAME FOR INSERTIONS
### A little trickier as for each length you need to pick the most common base at each insertion position (rather than a deletion where the base can only be reference...
	$min_indel_seq = "";
	@tempmaxbaseins = ();
	foreach $ins_length (sort {$a <=> $b} keys %{$minindel{"+"}}){
	$minindelfreq = $minindel{"+"}{$ins_length} / $depth;
	$minindelfreqrounded = sprintf("%.3f", $minindelfreq);
	$tempmaxins = 0;
	@tempmaxbaseins = ();
		$tempmaxinsbase = (reverse sort {$minindelseq{"+"}{$ins_length}{$a} <=> $minindelseq{"+"}{$ins_length}{$b} } keys %{$minindelseq{"+"}{$ins_length}})[0];
		$tempmaxins = $minindelseq{"+"}{$ins_length}{$tempmaxinsbase};
		foreach $base (keys %{$minindelseq{"+"}{$ins_length}}){
			if($minindelseq{"+"}{$ins_length}{$base} >= $tempmaxins){
			push(@tempmaxbaseins,$base);
			}
		}
	# if there are 2 or more bases with equal depths pick one at random.
	$tmbi = int rand (scalar @tempmaxbaseins);
	$min_indel_seq .= $tempmaxbaseins[$tmbi];
	
	$var = $ins_length.$min_indel_seq;
		if($minindelfreqrounded >= $mv_freq_cutoff){
			if($minindel{"+"}{$ins_length} >= $mv_variant_depth_cutoff){
			$identityrounded2 = sprintf("%.3f", ($identical/$depth));
			print MV "$seq\t$pos\t$refbase\t$depth\t".$majerror."MinIns\t$var\t".$minindel{"+"}{$ins_length}."\t$minindelfreqrounded\t0\t$fr_ratio_refbase\n";
			}
			if($minindelfreq >= 0.5){
				if($minindel{"+"}{$ins_length} < $mv_variant_depth_cutoff){
				print MV "$seq\t$pos\t$refbase\t$depth\tMinInsLow\t$var\t".$minindel{"+"}{$ins_length}."\t$minindelfreqrounded\t0\t$fr_ratio_refbase\n";
				$identityrounded2 = sprintf("%.3f", ($identical/$depth));
				print MV "$seq\t$pos\t$var\t$depth\tMajErrorMinInsLow\t$refbase\t$identical\t$identityrounded2\t0\t$fr_ratio_refbase\n";
				}
				else{
				$identityrounded2 = sprintf("%.3f", ($identical/$depth));
				print MV "$seq\t$pos\t$var\t$depth\tMajErrorMinIns\t$refbase\t$identical\t$identityrounded2\t0\t$fr_ratio_refbase\n";
				}
			$consensus_ins{$pos} = $var;
			}
		}
	}

	
	if($multidel != 0){
	$freq = $multidel / $depth;
	$rounded = sprintf("%.3f", $freq);
		if($rounded >= $mv_freq_cutoff){
			if($multidel >= $mv_variant_depth_cutoff){
			print MV "$seq\t$pos\t$refbase\t$depth\t".$majerror."MultiDel\t\*\t$multidel\t$rounded\t0\t$fr_ratio_refbase\n";
			}
			if($freq >= 0.5){
		#	print MV out new majority base in refbase position and refbase position in minority variant position plus the number of times reference base seen (i.e. .'s and ,'s summed)
		#	need this line to store new min var stats as %incorrect ref != 1-%majority variant necessarily does it...
		#	print LOG "md id $identical, $rounded, $rounded2\n";
				if($multidel < $mv_variant_depth_cutoff){
				print MV "$seq\t$pos\t$refbase\t$depth\tMultiDelLow\t\*\t$multidel\t$rounded\t0\t$fr_ratio_refbase\n";
			### need this if in cases where overall depth and variant depth are low but consensus is wrong.
			### e.g. depth 95, majority variant=48, min var=47
			### as the variant won't have printed out above...
				$rounded2 = sprintf("%.3f", ($identical/$depth));
				print MV "$seq\t$pos\t\*\t$depth\tMajErrorMultiDelLow\t$refbase\t$identical\t$rounded2\t0\t$fr_ratio_refbase\n";
				}
				else{
				$rounded2 = sprintf("%.3f", ($identical/$depth));
				print MV "$seq\t$pos\t\*\t$depth\tMajErrorMultiDel\t$refbase\t$identical\t$rounded2\t0\t$fr_ratio_refbase\n";
				}
			$multidel_pos_check{$pos}++;
			}
		}
	$acgtd{$pos}{"MD"} = $multidel;
	}

}#end if  depth>depthcut statement
}


####testing####

###$consensus_del{"500"} = "10Acgtacgtag";



foreach $pos (sort {$a<=>$b} keys %consensus_del){
$var = $consensus_del{$pos};
	if($var =~ /(\d+)(\S+)/){
	$x = 0;
	$varlen = $1;
	$bases = $2;
		while($x < $varlen){
		$consensus_sequence{($pos+$x+1)} = "";
		$x++;
		}
	}
}


###HAVE TO DO THIS - MPILEUP SCREWS UP WHEN DELETIONS BETWEEN CIGAR PADS APPARENTLY - DON'T GET THE -1A, -1T etc IN THE PREVIOUS PILEUP LINE BUT STILL GET 8000 *'s FOR EXAMPLE... UUGHGHGHGHGH
foreach $pos (sort {$a<=>$b} keys %multidel_pos_check){
$consensus_sequence{$pos} = "";
}



###testing
###$consensus_ins{"515"} = "10dgfdgfdgfd";

foreach $pos (sort {$a<=>$b} keys %consensus_ins){
#nope do this in output maybe 
$var = $consensus_ins{$pos};
	if($var =~ /(\d+)(\S+)/){
	$x = 0;
	$varlen = $1;
	$bases = $2;
		while($x < $varlen){
		$consensus_sequence{($pos+(0.0001*($x+1)))} = substr($bases,$x,1);
		$consensus_depth{($pos+(0.0001*($x+1)))} = $consensus_depth{$pos};
		$x++;
		}
	}
}


## Working out sliding window depth value for each position
## Am including repeat values for insertions as if there's a long insertion with v high depth it seems
## daft to penalise downstream bases where the original gapfilled sequence say is wrong with 0 coverage
## but there's a 40bp insertion
## to make amends for it...

@allpos = (sort {$a<=>$b} keys %consensus_depth);
$i = 0;
$maxpos = $#allpos;

foreach $pos (sort {$a<=>$b} keys %consensus_depth){
@window = ();
@sorted_window = ();
$bottom = $i - int ($sliding_window_size/2);
$top = $i + int ($sliding_window_size/2);
	if($bottom < 0){
	$bottom = 0;
	}
	if($top > $maxpos){
	$top = $maxpos;
	}
	for($j = $bottom; $j <= $top ; $j++){
	push(@window,$consensus_depth{$allpos[$j]});
	}
@sorted_window = sort {$a <=> $b} (@window);
	if($#sorted_window/2-int($#sorted_window/2) == 0){
	$window_depth{$pos} = $sorted_window[int($#sorted_window/2)];
	}
	##almost unnecessary as script stands as there will always be an odd number of depths in window...
	else{
	$window_depth{$pos} = 0.5 * ($sorted_window[($#sorted_window/2)] + $sorted_window[(($#sorted_window/2)-1)]);
	}
$i++;
}




open(LOG_CONS,">$logcons") || die $!;
$final_consensus = "";
foreach $pos (sort {$a<=>$b} keys %consensus_sequence){
print LOG_CONS "$pos\t".$consensus_sequence{$pos}."\t".$consensus_depth{(int $pos)}."\t".$window_depth{$pos}."\n";
	#if($consensus_depth{$pos} > $cons_depth_cutoff){
	#doing either or as just using median messes up good edge bases when a more internal base has 99 or less...
	
	if($consensus_depth{$pos}>0){
		if(($consensus_depth{$pos} > $cons_depth_cutoff) || ($window_depth{$pos} > $cons_depth_cutoff)){
		$final_consensus .= $consensus_sequence{$pos};
		}
# see if we can eke out sequence for low coverage parts - eg depth=5 normally => N, but smalt is mapping 4 reads across identically
# so perhaps we should keep the depth of 4 consensus part rather than N's. 
# alternative is to alter majvarerror2.pl to grep -v Low before grep -c Maj bit...

		else{
		$final_consensus .= lc($consensus_sequence{$pos});
		}
	}
	else{
	$final_consensus .= "N";
	}
}
close(LOG_CONS) || die $!;

### This is from pre-iteration era - with each iteration the word "consensus" gets added every time = daft!
###$consensus_genome_name = $refname.".consensus";
if($refname =~ /(\S+)\.genome/){
$consensus_genome_name = $1.".consensus1";
}
elsif($refname =~ /(\S+)\.consensus(\d+)/){
$consensus_genome_name = $1.".consensus".($2+1);
}
elsif($refname =~ /(\S+)/){
$consensus_genome_name = $1.".consensus";
}

print CONS ">$consensus_genome_name\n";
print CONS "$final_consensus\n";


close(IN) || die $!;
close(CONS) || die $!;
close(MV) || die $!;


## NB Altered loop to take care of majority insertions - currently the position gets upped by 0.0001
## to make outputting the consensus a doddle (sorting by position). HOWEVER - this screws up the %acgtd keys
## in the loop below - need to change 10, 11, 11.0001, 11.0002, 12 to 10,11,12,13,14 etc...
## TODO: tidy this up by adding inserted bases and 0.0001ish positions into acgtd in consensus_ins 
## loop line 575ish above.  Still would need to do the pos - int pos thing below so not ideal either way.

open(ACGTD,">$base_freq_out") || die $!;
print ACGTD "Position\tA\tC\tG\tT\tDel\n";
$posadj = 0;
@columns = qw(A C G T MD);
	foreach $pos (sort {$a <=> $b} keys %consensus_depth){
		if(($pos - (int $pos)) != 0){
		$base = $consensus_sequence{$pos};
        	$acgtd{$pos}{$base} = $consensus_depth{$pos};
		$posadj++;
        	}
        $acgtdpos=(int $pos)+$posadj;
	
	print ACGTD $acgtdpos;
		foreach $column (@columns){
			if(!exists($acgtd{$pos}{$column})){
			print ACGTD "\t0";
			}
			else{
			print ACGTD "\t".$acgtd{$pos}{$column};
			}
		}
	print ACGTD "\n";
	}

close(ACGTD) || die $!;


