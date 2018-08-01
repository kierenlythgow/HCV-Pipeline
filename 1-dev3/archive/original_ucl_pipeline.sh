#!/bin/bash -l
#$ -S /bin/bash
#$ -l h_rt=12:00:0
#$ -l tmpfs=25G
#$ -l mem=2G
#$ -pe smp 8
#$ -N name_of_job_for_cluster
#$ -wd /path/to/working_directory
#$ -t 1-400


# For a full run you only need to specify the results folder and in a tab-delimited file
# a list of sample IDs and the path to their forward and reverse reads 
# (these can be gzipped if necessary).

ResultsFolder=/home/kieren/UCL_IVA ; 
#hiv
paramfile=/home/kieren/UCL_IVA/hiv_samples.txt
#hev
#paramfile=/home/kieren/UCL_IVA/hev_samples.txt
# e.g. of paramfile: tab-delimited text file
# sample1	/path/to/sample1_fwd_reads.fq	/path/to/sample1_rev_reads.fq
# sample2	/path/to/sample2_fwd_reads.fq	/path/to/sample2_rev_reads.fq
# etc.

#Added 17/11/2016 by KL
#number=$1


sample=`sed -n ${number}p $paramfile | awk '{print $1}'`
Reads1=`sed -n ${number}p $paramfile | awk '{print $2}'`
Reads2=`sed -n ${number}p $paramfile | awk '{print $3}'`


# Set the working directory to either $TMPDIR on the cluster or $wd on desktop

##################################
## For cluster                  ##
## i.e. commment out on desktop ##
#   cd $TMPDIR ; 
# # # OR # # # # # # # # #  # # ##
## For desktop                  ##
## i.e. comment out on Legion   ##
cd /home/kieren/UCL_IVA        #                
##################################


#Original
#mkdir -p output/temp;
#NB IVA fails if /path/to/working_directory/output/assembly exists
#mkdir logs;

#mkdir -p $ResultsFolder/$sample/output/temp ; 
#mkdir -p $ResultsFolder/$sample/output/assembly ;
#mkdir -p $ResultsFolder/$sample/logs ;
#cd output;

#KL changes
mkdir -p $ResultsFolder/$sample/output/temp ;
mkdir -p $ResultsFolder/$sample/output/assembly ;
mkdir -p $ResultsFolder/$sample/logs ;
cd $ResultsFolder/$sample/output;

###TRIM READS###
java -jar /home/kieren/UCL_IVA/trimmomatic-0.33.jar PE -threads 8 $Reads1 $Reads2 "temp/"$sample"_hq_1.fastq" "temp/"$sample"_unpaired_1.fastq" "temp/"$sample"_hq_2.fastq" "temp/"$sample"_unpaired_2.fastq" ILLUMINACLIP:/path/to/adapters.fasta:2:10:7 LEADING:10 TRAILING:10 SLIDINGWINDOW:4:30 MINLEN:50 ;
#rm temp/*unpaired_*.fastq ; 
#cp "temp/"*hq*.fastq $ResultsFolder/$sample"/output/temp/". ; 
#rm temp/*hq_*.fastq ; 

###DECOY STEP START###
#hiv
smalt map -x -y 0.5 -i 500 -n 8 /home/kieren/UCL_IVA/hg38_hiv_k15_s3 "temp/"$sample"_hq_1.fastq" "temp/"$sample"_hq_2.fastq" | awk '{if ($3 !~ /^chr|\*/ && $7 !~ /^chr|\*/) print $0}' > "temp/"$sample"_pairs.sam" ;
#hev
#smalt map -x -y 0.5 -i 500 -n 8 /home/kieren/UCL_IVA/hg38_hev_k15_s3 "temp/"$sample"_hq_1.fastq" "temp/"$sample"_hq_2.fastq" | awk '{if ($3 !~ /^chr|\*/ && $7 !~ /^chr|\*/) print $0}' > "temp/"$sample"_pairs.sam" ;
samtools view -bhf 64 "temp/"$sample"_pairs.sam" | samtools bam2fq -  > $sample"_filtered_1.fastq" ;
samtools view -bhf 128 "temp/"$sample"_pairs.sam" | samtools bam2fq -  > $sample"_filtered_2.fastq" ;
#rm temp/*.sam;
gzip -c $sample"_filtered_1.fastq" > $sample"_filtered_1.fastq.gz" ; 
gzip -c $sample"_filtered_2.fastq" > $sample"_filtered_2.fastq.gz" ; 
### Copy filtered reads across to results folder in case job takes ages (e.g. HiSeq) as otherwise lose when wall clock time expires...###
cp $sample"_filtered_1.fastq.gz" $ResultsFolder/$sample"/output/"$sample"_filtered_1.fastq.gz" ;
cp $sample"_filtered_2.fastq.gz" $ResultsFolder/$sample"/output/"$sample"_filtered_2.fastq.gz" ;

###Starting IVA###
#For now, leave out trimmomatic from IVA - time is of the essence. Can always repeat and see if it makes a difference...
iva --max_contigs 50 -t 8 -f $sample"_filtered_1.fastq" -r $sample"_filtered_2.fastq" assembly   ;

# OR if you have primers: I didn't for this example
# iva --trimmomatic /path/to/trimmomatic-0.33.jar -pcr_primers /path/to/primers.fasta --adapters /path/to/adapters.fasta --max_contigs 50 -t 8 -f $sample"_filtered_1.fastq" -r $sample"_filtered_2.fastq" assembly   ;
# # # # 


#Assuming IVA works, rename IVA contigs to sample specific contigs
sed -r {'s/>contig/>'$sample'/'} assembly/contigs.fasta > "assembly/"$sample".contigs.fasta" ; 
## Copy contigs across to results folder in case job has taken ages - otherwise lose all data when wall clock time is up...###
cp "assembly/"$sample".contigs.fasta" $ResultsFolder/$sample"/output/assembly"/. ;
### IVA done - now finding best reference###
#hiv
/home/kieren/lastz-distrib/bin/lastz "assembly/"$sample".contigs.fasta"[multiple] /home/kieren/UCL_IVA/Data/hiv1.fasta --ambiguous=iupac --format=GENERAL > "assembly/"$sample".contigs.lastz" ; 
#hev
#/home/kieren/lastz-distrib/bin/lastz "assembly/"$sample".contigs.fasta"[multiple] /home/kieren/UCL_IVA/Data/HEV_orf2_BBVUprimers.fasta --ambiguous=iupac --format=GENERAL > "assembly/"$sample".contigs.lastz" ;
#hiv
perl -s /home/kieren/UCL_IVA/Command_line_pipeline/lastz_bestref.pl -contig_lastz="assembly/"$sample.contigs.lastz  -blastdb=/home/kieren/UCL_IVA/Data/hiv1.fasta -best_ref_fasta=$sample"-ref.fasta" -lastz_best_hit_log="../logs/"$sample"_lastz_besthit.log" ; 
#hev
#perl -s /home/kieren/UCL_IVA/Command_line_pipeline/lastz_bestref.pl -contig_lastz="assembly/"$sample.contigs.lastz  -blastdb=/home/kieren/UCL_IVA/Data/HEV -best_ref_fasta=$sample"-ref.fasta" -lastz_best_hit_log="../logs/"$sample"_lastz_besthit.log" ;

### Have best reference, now assemble draft genome ###
/home/kieren/lastz-distrib/bin/lastz $sample"-ref.fasta" "assembly/"$sample".contigs.fasta" --ambiguous=iupac > "assembly/"$sample".contigs-vs-bestref.lav" ; 
perl -w -s /home/kieren/UCL_IVA/Command_line_pipeline/lastz_analyser.WITH_REVCOMP.pl -reference_fasta_file=$sample"-ref.fasta" -sample_fasta_file=assembly/$sample".contigs.fasta" -lastz_results_file="assembly/"$sample".contigs-vs-bestref.lav" -cutoff=50000 -with_revcomp=yes -output=temp/$sample".lastz_analysed_file" -log_file="../logs/"$sample"_lastz_analyser.log" ;
## map reads to contigs to work out which sequences to choose when contigs overlap
smalt index -k 15 -s 3 "assembly/"$sample".contigs.k15_s3" "assembly/"$sample".contigs.fasta" ; 
smalt map -x -y 0.5 -i 500 -n 8 -f bam -o "temp/"$sample".contigs.bam" "assembly/"$sample".contigs.k15_s3" $sample"_filtered_1.fastq" $sample"_filtered_2.fastq" ; 
samtools sort -@ 8 -T temp_sort -o "temp/"$sample".contigs.sorted.bam"  "temp/"$sample".contigs.bam"
samtools index "temp/"$sample".contigs.sorted.bam"
samtools mpileup -f "assembly/"$sample".contigs.fasta" -d 1000000 -o "temp/"$sample".contigs.mpileup" "temp/"$sample".contigs.sorted.bam";
awk {'print $1"\t"$2"\t"$3"\t"$4'} "temp/"$sample".contigs.mpileup" > $sample".contigs.fake_mpileup" ; 
perl -w -s /home/kieren/UCL_IVA/Command_line_pipeline/genome_maker2b.pl -sample_pileup_file="temp/"$sample".contigs.mpileup" -contigs="assembly/"$sample".contigs.fasta" -reference_mapped_consensus=$sample"-ref.fasta" -lastz_analysed_file="temp/"$sample".lastz_analysed_file" -ref_correct_start=0 -ref_correct_stop=20000 -output="temp/"$sample"-genome.fasta" -logfile="../logs/"$sample"_genome_maker.log" ;
#rm temp/*.bam ;
#rm temp/*.mpileup ;
##now we've made the draft genome, time to map all the reads onto it and make the first consensus sequence...	
smalt index -k 15 -s 3 "temp/"$sample"-genome.k15_s3" "temp/"$sample"-genome.fasta" ; 
smalt map -x -y 0.5 -i 500 -n 8 -f bam -o "temp/"$sample"-genome.bam" "temp/"$sample"-genome.k15_s3" $sample"_filtered_1.fastq" $sample"_filtered_2.fastq" ; 
samtools sort -@ 8 -T temp_sort -o "temp/"$sample"-genome.sorted.bam"  "temp/"$sample"-genome.bam"
samtools index "temp/"$sample"-genome.sorted.bam"
samtools mpileup -f "temp/"$sample"-genome.fasta" -d 1000000 -o "temp/"$sample"-genome.mpileup" "temp/"$sample"-genome.sorted.bam";
awk {'print $1"\t"$2"\t"$3"\t"$4'} "temp/"$sample"-genome.mpileup" > "temp/"$sample"-genome.fake_mpileup" ;
#rm temp/*.bam;
perl -w -s /home/kieren/UCL_IVA/Command_line_pipeline/cons_mv.pl -mpileup="temp/"$sample"-genome.mpileup" -reference_fasta="temp/"$sample"-genome.fasta" -mv_freq_cutoff=0.01 -mv_overall_depth_cutoff=100 -mv_variant_depth_cutoff=20 -cons_depth_cutoff=80 -sliding_window_size=300 -consensus_out="temp/"$sample".consensus1.preNcut.fasta" -mv_out="temp/"$sample"-genome.fasta.mv" -base_freq_out="temp/"$sample"-genome.fasta.basefreqs.tsv" ;
perl -w -s /home/kieren/UCL_IVA/Command_line_pipeline/N_remover_from_consensus.pl -cutoff=46 "temp/"$sample".consensus1.preNcut.fasta" > "temp/"$sample".consensus1.fasta";
#rm temp/*.mpileup;
##experience shows that one round of iteration is rarely enough so do one more...
smalt index -k 15 -s 3 "temp/"$sample".consensus1.k15_s3" "temp/"$sample".consensus1.fasta" ; 
smalt map -x -y 0.5 -i 500 -n 8 -f bam -o "temp/"$sample".consensus1.bam" "temp/"$sample".consensus1.k15_s3" $sample"_filtered_1.fastq" $sample"_filtered_2.fastq" ;
samtools sort -@ 8 -T temp_sort -o "temp/"$sample".consensus1.sorted.bam"  "temp/"$sample".consensus1.bam" ; 
samtools index "temp/"$sample".consensus1.sorted.bam" ; 
samtools mpileup -f "temp/"$sample".consensus1.fasta" -d 1000000 -o "temp/"$sample".consensus1.mpileup" "temp/"$sample".consensus1.sorted.bam";
awk {'print $1"\t"$2"\t"$3"\t"$4'} "temp/"$sample".consensus1.mpileup" > "temp/"$sample".consensus1.fake_mpileup" ;
#rm temp/*.bam;
perl -w -s /home/kieren/UCL_IVA/Command_line_pipeline/cons_mv.pl -mpileup="temp/"$sample".consensus1.mpileup" -reference_fasta="temp/"$sample".consensus1.fasta" -mv_freq_cutoff=0.01 -mv_overall_depth_cutoff=100 -mv_variant_depth_cutoff=20 -cons_depth_cutoff=80 -sliding_window_size=300 -consensus_out="temp/"$sample".consensus2.preNcut.fasta" -mv_out="temp/"$sample".consensus1.fasta.mv" -base_freq_out="temp/"$sample".consensus1.fasta.basefreqs.tsv" ;
perl -w -s /home/kieren/UCL_IVA/Command_line_pipeline/N_remover_from_consensus.pl -cutoff=46 "temp/"$sample".consensus2.preNcut.fasta" > "temp/"$sample".consensus2.fasta";
#rm temp/*.mpileup;		
perl -w -s /home/kieren/UCL_IVA/Command_line_pipeline/majvarcheck2.pl -mvpath="temp/"$sample".consensus1.fasta.mv" -basefreq="temp/"$sample".consensus1.fasta.basefreqs.tsv" -fwdreads=$sample"_filtered_1.fastq" -revreads=$sample"_filtered_2.fastq";	

### All done - now copy results from working directory to results directory
### This is essential when running a job on the cluster, otherwise can do at your
### own leisure...
cp temp/*.consensus.* . ;
cp temp/*.consensus.* $ResultsFolder/$sample/output/. ;
cd .. ;
cp -pr logs/* $ResultsFolder/$sample/logs/. ;
cp -pr output/temp/* $ResultsFolder/$sample/output/temp/. ;
cp -pr output/assembly/* $ResultsFolder/$sample/output/assembly/. ; 
