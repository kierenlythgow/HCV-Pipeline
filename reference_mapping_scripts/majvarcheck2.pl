#!/usr/bin/perl -s

## Start iteration from second consensus sequence which should be named 12345_1_1.consensus2.fasta and the minority variants file 12345_1_1.consensus1.fasta.mv
## -mvpath=OUTPUT_PATH/14535_1_12.consensus1.fasta.mv
## -fwdreads=PATH_TO_QC_READS/sampleQCreads.f.fq
## -revreads=PATH_TO_QC_READS/sampleQCreads.r.fq
## -basefreq=path to the base_freq_out file with ACGTDel counts for each position (created by cons_mv.pl)


#$mv = $ARGV[0];

use File::Basename;
use Cwd 'abs_path';

if(!$fwdreads){
die "No FASTQ file for fwd reads given\n";
}

if(!$revreads){
die "No FASTQ file for rev reads given\n";
}

if(!$mvpath){
die "Error: No mvpath given\n";
}

if(!$basefreq){
die "Error: No base frequencies file given\n";
}
else{
$base_freq_out = $basefreq;
}

$majvarcount = 10000;
$majvarcount_OLD = 10000;
$iterations = 0;
while((($majvarcount+$majvarcount_OLD) > 0) && ($iterations<10)){
  $majvarcount = `grep -c Maj $mvpath`;
  chomp($majvarcount);
	if($mvpath =~ /(\S+sus)(\d+)\.fasta/){
  		if($2>1){
		$mvpath_OLD = $1.($2-1).".fasta.mv";
		}
		else{
		$mvpath_OLD = $mvpath;
		$mvpath_OLD =~ s/\.consensus1/\-genome/;
		}
	}
#print "DAN:mvpath:$mvpath\tmvpath_OLD:$mvpath_OLD\n";
  $majvarcount_OLD = `grep -c Maj $mvpath_OLD`;
  chomp($majvarcount_OLD);
#print "DAN:majvarcount:$majvarcount\tmajvarcount_OLD:$majvarcount_OLD\n\n";
  if(($majvarcount+$majvarcount_OLD) > 0){
    if($mvpath =~ /(\S+sus)(\d+)\.fasta/){
      $x = $2;
      $x++;
      $mvnew = $1.$x.".fasta.mv";
      $y = $x + 1;
      $smalt_index = $1.$x.".fasta.index";
      $consensus_fasta = $1.$x.".fasta";
      $sam = $1.$x.".fasta.smalt.sam";
      $bam = $1.$x.".fasta.smalt.bam";
      $sorted = $1.$x.".fasta.smalt.sorted";
      $sorted_bam = $1.$x.".fasta.smalt.sorted.bam";
      $pileup = $1.$x.".fasta.smalt.sorted.pileup";
      $preNcut_fasta = $1.$y.".preNcut.fasta";
      $new_consensus_fasta = $1.$y.".fasta";
      $base_freq_out = $1.$x.".fasta.basefreqs.tsv";

      $final_consensus = $1.".fasta";
      $final_mv = $1.".mv";
      $final_base_freq = $1.".basefreqs.tsv";
      $final_fake_mpileup = $1.".fake_mpileup";

      system("smalt index -k 15 -s 3 $smalt_index $consensus_fasta");
      system("smalt map -x -y 0.5 -i 500 -n 8 -o $sam $smalt_index $fwdreads $revreads");
      system("samtools view -@ 8 -bh -S $sam > $bam");
      system("samtools sort -@ 8  $bam $sorted");
      system("samtools faidx $consensus_fasta");
      system("samtools mpileup -f $consensus_fasta $sorted_bam -d 1000000 > $pileup");
	$fake_mpileup = $pileup;
	$fake_mpileup =~ s/pileup$/fake_mpileup/;
	$fake_mpileup =~ s/\.smalt\.sorted//;
      system('awk {\'print $1"\t"$2"\t"$3"\t"$4\'} '."$pileup > $fake_mpileup");
	$dirname  = dirname(abs_path($0));
      system("perl -w -s $dirname/cons_mv.pl -mpileup=$pileup -reference_fasta=$consensus_fasta -consensus_out=$preNcut_fasta -mv_out=$mvnew -base_freq_out=$base_freq_out");
      system("perl -w -s $dirname/N_remover_from_consensus.pl -cutoff=46 $preNcut_fasta > $new_consensus_fasta");
      system("rm $sam $bam $sorted_bam $pileup");
      }
    $mvpath = $mvnew;
    $iterations++;
  }
}

### Now tidy up - renaming final iteration as .consensus.fasta, renaming final pileup and bam file  and then removing all other pileup, sam, bam, thank-you-mam, etc.

if($new_consensus_fasta){
  system("cp $new_consensus_fasta $final_consensus");
  system("cp $mvpath $final_mv");
  system("cp $base_freq_out $final_base_freq");
  system("cp $fake_mpileup $final_fake_mpileup");
  system("sed -r {'s/consensus.*/consensus/'} $new_consensus_fasta > $final_consensus") ;
}
else{
  if($mvpath =~ /(\S+consensus\d+\.fasta)\.mv/){
    $first_consensus = $1;
    $final_consensus = $first_consensus;
    $final_consensus =~ s/sus\d+/sus/;
    system("cp $first_consensus $final_consensus");
    $final_mv = $mvpath;
    $final_mv =~ s/sus\d+\.fasta/sus/;
    $final_base_freq = $base_freq_out;
    $final_base_freq =~ s/sus\d+\.fasta/sus/;
    $fake_mpileup = $first_consensus;
    $fake_mpileup =~ s/\.fasta$/\.fake_mpileup/;
    $final_fake_mpileup = $final_consensus;
    $final_fake_mpileup =~ s/sus\.fasta/sus\.fake_mpileup/;
    system("cp $fake_mpileup $final_fake_mpileup");
    system("cp $mvpath $final_mv");
    system("cp $base_freq_out $final_base_freq");
    system("sed -r {'s/consensus.*/consensus/'} $first_consensus > $final_consensus") ;
  }
}
