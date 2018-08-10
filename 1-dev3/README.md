###README
======

The generate_hcv_consensus component performs a series of processes that consumes trimmed fastq files and generates a consensus sequence and minority variant report.

The processes involve:

1. Human filtering - removal of human reads from the sample.
2. Split population analysis - this aims to detect the presence of mixed infections.
3. De novo assembly - VICUNA is used to generate contigs.
4. Optimum reference selection - an optimum reference sequence is selected from a database of HCV sequences.
5. Draft assembly - Using the optimal reference sequence a draft assembly is created.
6. Contig mapping - contigs are mapped against the draft assembly.
7. Consensus & minority variants - Quasibam generates a consensus and minority variant report.
8. Combine XMLS - XML's are combined into one document.

###REQUIREMENTS

  1) python/2.7.6                                                   
  2) yaml/1.1                                                      
  3) psutil/python2.7.6/3.0.0                                      
  4) jdk/1.8.0_121                                                 
  5) picard-tools/1.111                                            
  6) samtools/1.1                                                  
  7) blast+/2.2.27                                                 
  8) bamtools/2.3.0
  9) bwa/0.7.13
  10) smalt/0.7.6
  11) bedtools/2.26.0
  12) vicuna/1.3
  13) phe/common_modules/1-22
  14) phe/quasi_bam

###PREREQUISITES

generate_hcv_consensus.py -i <input_dir> -w <workflow> -r <reference_sets>

###HISTORY

Contributors

Kieren Lythgow, Matthew Goulden
