**README**
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

**REQUIREMENTS**

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

**PREREQUISITES**

```generate_hcv_consensus.py -i <input_dir> -w <workflow> -r <reference_sets>```

**RETURN CODES**

    66: no input provided in command line 
    67: component directorty already exists
    68: no processed fastqs in the input directory
    69: failed to unzip .gz fastq files
    70: paired sam file does not exist
    71: paired sam file is empty 
    72: samtools view failed 
    73: samtools bam2fq failed 
    74: fastq2sam failed 
    75: filtered bam file does not exist 
    76: filtered bam file is empty 
    77: snork splitpops failed
    78: GENERATE_HCV_CONSENSUS_VICUNA_CONFIG environment variable not set
    79: vicuna de novo assembly failed
    80: lastz failed
    81: lastz bestref failed
    82: lastz compare failed
    83: lastz analyser reverse failed
    84: bwa index failed
    85: bwa mem failed
    86: samtools view failed in contig mapping
    87: samtools sort failed
    88: samtools index failed
    89: samtools mpileup failed
    90: genome maker failed
    91: cons mv failed
    92: n remover failed
    93: cons mv basefreq failed
    94: n remover cons2 failed
    95: majvar bwa failed
    96: samtools view final failed
    97: quasibam failed
    98: expecting 4 xml files


**HISTORY**

**Contributors**

Kieren Lythgow, Matthew Goulden
