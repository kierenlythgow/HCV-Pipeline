#!/usr/bin/perl -s

if(!$contig_lastz){
die "No contig lastz input file given: use -contig_lastz=PATH/PILEUP_FILE\nMake using e.g. lastz sample.all_contigs.fasta[multiple] hiv-1.lastz.fasta[multiple] --ambiguous=IUPAC --format=GENERAL > contigs.lastz.file\n\n";
}
if(!$blastdb){
$blastdb = '/path/to/BLAST/virus/reference/database/BLAST_viral_DB';
}
if(!$best_ref_fasta){
$best_ref_fasta = $contig_lastz;
$best_ref_fasta =~ s/\..*/\-ref\.fasta/;
	if($best_ref_fasta eq $contig_lastz){
	$best_ref_fasta = $contig_lastz.".ref.fasta";
	}
}

if(!$lastz_best_hit_log){
$lastz_best_hit_log = "logfile.".$best_ref_fasta.".log";
}
###
###
# REMEMBER - Because Baz didn't trim the non-coding ends off each reference sequence,
# a full length subtype D sample say (ie more than just CWG) can give a better
# BLAST alignment score to full length hxb2 than a CWG subtype D reference
# OUCH!!!!!!!!!!!!
# That's not a big problem, but in instances where sequencing is gappy
# but covers the 5' LTR say, this could be a big issue for gap-filling.
# e.g. that comparison where we used diff subtypes and got diff envelope...
# So in future we need to use two databases I think, one CWG for the lastz/BLAST step
# and then look up the full sequence for the rest of the pipeline...
###
###
#
#
# NB loads of commenting out - have decided to adapt a partially finished script
# to mimic Baz's BLASTN fudge - he just picks the top score
# so for consistency that's all I'll do.
###

$fudge{"0"}{"NoLastzHit"} = 1;
open(CONTIG_LASTZ,"$contig_lastz") || die $!;
while(<CONTIG_LASTZ>){
##e.g.
##score  name1   strand1 size1   zstart1 end1    name2   strand2 size2   zstart2 end2    identity        idPct   coverage        covPct
#6669    15228_1_84.1    +       3291    248     324     B#A04321.1      +       9193    9117    9193    72/76   94.7%   76/3291 2.3%
#269257  15228_1_84.1    +       3291    248     3291    B#A04321.1      +       9193    21      3064    2918/3041       96.0%   3043/3291 92.5%

	if(/^(\d+)\s+.*\n/){
	$score = $1;
	chomp($_);
	@line = split(/\t/,$_);
	$contig = $line[1];
	$contig_length = $line[3];
	$ref = $line[6];
	$reflength = $line[8];
	#need to decide whether to pick best ref based on how much of contig is covered or by longest ref alignment...
	#think it should probably be based primarily on score, then by length of ref in alignment, then by seqID of ref in alignment
	$refstart=$line[9];
	$refstop=$line[10];
	$ids_of_ref = $line[11];
	$id = $line[12];
	$cov_of_contig = $line[13];
	$covpct = $line[14];
	#will be picking by top score
	$score{$contig}{$score}{$ref} = $_;
	#store all hits so can search for other alignments to same ref later on e.g. for gappy contigs 
	$all_hits{$contig}{$ref}{$_}++;
	# this will take ages for 5459 hits x 50 contigs... go with best score for now
	#	for($i=$refstart; $i<=$refstop; $i++){
	#	$refcov{$ref}{$contig}{$i}++;
	#	}
	$fudge{$score}{$ref}++;
	}
}
close(CONTIG_LASTZ) || die $!;
open(LOG,">$lastz_best_hit_log") || die $!;
@scores = reverse sort {$a <=> $b} (keys %fudge);
$top_score=$scores[0];
@best_refs = keys %{$fudge{$top_score}};
	if(scalar @best_refs > 1){
	$index = int(rand(scalar @best_refs));
	$best_refname = $best_refs[$index];
	print LOG "MultipleTopHits:",@best_refs,"\n";
	}
	else{
	$best_refname = $best_refs[0];
	}
$blastdbcmd = "blastdbcmd -db $blastdb -dbtype nucl -entry \"$best_refname\" > $best_ref_fasta";
print LOG "TopHit:$best_refname\nTopScore:$top_score\n";
print LOG $blastdbcmd,"\n";
	if($best_refname ne "NoLastzHit"){
	system($blastdbcmd);
	}
	else{
	print LOG "ERROR: lastz was unable to find a single contig alignment to references\n\n";
	}
close(LOG) || die $!;

#open(CONTIG_PILEUP,"$contig_pileup") || die $!;
#while(<CONTIG_PILEUP>){
##i.e. 12559_1_10#1    477     G       6042
#	if(/^(\S+)\s+(\d+)\s+\S+\s(\d+)/){
#	$contig = $1;
#	$pos = $2;
#	$depth = $3;
#	$depth_hash{$contig}{$pos}=$depth;
#	}
#
#}
#close(CONTIG_PILEUP) || die $!;

#foreach $contig (keys %depth_hash){
#	foreach $pos (keys %{$depth_hash{$contig}){
#	$depthsum{$contig}+=$depth_hash{$contig}{$pos};
#	$poscount{$contig}++;
#	}
#}

#foreach $contig (sort keys %score){
#@temp_scores = reverse sort {$a <=> $b} (keys %score);
#$top_score = $temp_scores[0];
#@temp_refs = sort keys %{$score{$contig}{$top_score}};
#	if(scalar @temp_refs > 1){
#	#Houston, we have a problem of tied references
#	%refalignlength = ();
#		foreach $ref (@temp_refs){
#		$hitline = $score{$contig}{$score}{$ref};
#		@hitline = split(/\t/,$hitline);
#		$temp_refalignlength = $hitline[10]-$hitline[9]+1;
#		$refalignlength{$temp_refalignlength}{$ref}++;
#		}
#	@align_lengths = reverse sort {$a<=>$b} (keys %refalignlength);
#	$top_align_length = $align_lengths[0];
#	@temp_refs = sort (keys %{$refalignlength{$top_align_length}});
#		if(scalar @temp_refs > 1){
#		##argh - we now have to go down to sequence identity
#			foreach $ref (@temp_refs){
#			$hitline = $score{$contig}{$score}{$ref};
#			@hitline = split(/\t/,$hitline);
#			$temp_seqIDref = $hitline[12];
#			$seqID{$temp_seqIDref}{$ref}++;
#			}
#		@seqID = reverse sort {$a <=> $b} (keys %seqID);
#		$top_seqID = $seqID[0];
#		@temp_refs = sort (keys %{$seqID{$top_seqID}});
#			if(scalar @temp_refs > 1){
#			$index = scalar rand(scalar @temp_refs);
#			$best_ref{$contig} = $temp_refs[$index];
#			}
#			else{
#			$best_ref{$contig} = $temp_refs[0];
#			}
#		}
#		else{
#		$best_ref{$contig} = $temp_refs[0];
#		}
#	}
#
#	else{
#	$best_ref{$contig}=$temp_refs[0];
#	}
#}