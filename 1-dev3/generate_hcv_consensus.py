#!/usr/bin/env python

"""Need to ensure a folder exists for the particular sample containing the trimmed fastqs.
This is required as VICUNA's input is a directory that needs to contain only fastqs that
require assembly."""

import os
import sys
import glob
import gzip
import subprocess
import argparse
import re
import logging
import csv
import shutil
from xml.dom import minidom
import xml.etree.ElementTree as ET
import inspect
import log_writer
from utility_functions import write_component_complete

__version__ = '0.1'
__date__ = '01Aug2018'
__author__ = 'kieren.lythgow@phe.gov.uk'

def parse_args():

    """
     Parge arguments
     Parameters
     ----------
     no inputs
     Returns
     -------
     Args: obj
     arguments object
     """

    Description = 'version %s, date %s, author %s' %(__version__, __date__, __author__)

    Parser = argparse.ArgumentParser(description=Description)

    Parser.add_argument('-i', '--input',
                          metavar='STRING',
                          dest='input',
                          required=True,
                          help='Input folder incl. path. [Required.]')

#     Parser.add_argument('-r', '--resultxml',
#                          metavar='STRING',
#                          dest='resxml',
#                          required=True,
#                          help='Name of result XML file. [Required.]')

    Parser.add_argument('-w', '--workflow',
                          metavar='STRING',
                          dest='workflow',
                          default='0',
                          help='The workflow that this sample was entered into. \
                                [default: 0 = get from file name]')

    Parser.add_argument('-ref', '--ref_set',
                        dest='refset_dir',
                        default = os.environ.get('GENERATE_HCV_CONSENSUS_REFERENCE_DIR'),
                        help='Full path to refset_dir including all reference requirements for this component [TODO succinct descr required]')

    Args = Parser.parse_args()
    return Args


def main():

    """
    Main function
    Parameters
    ----------
    no inputs
    Returns
    -------
    0 on success
    else int > 0
    """

    '''Return codes ------------ 
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
    '''


    Args = parse_args()
    if not Args.refset_dir:
        print('pass a full refset path to the `ref_set` passed to param')
        sys.exit(inspect.currentframe().f_lineno)

    #set up logger
    logger = log_writer.setup_logger('logs/generate_hcv_consensus.stdout',
                                     'logs/generate_hcv_consensus.stderr')
    log_writer.info_header(logger, "starting generate_hcv_consensus")
    log_writer.write_log(logger,
                         "processing samples in %s ...\n" % Args.input,
                         "info")

    #Check input is provided
    if not Args.input:
        log_writer.write_log(logger,'No input provided via --input/-i. '\
                'Nothing to do. Exiting.'
                'error')
        # electing to have a non-zero exit here although could be run like
        # this to find version from printed output i.e., not really an error
        sys.exit(66)

    #Determine sample name
    input_dir = Args.input
    
    #Change to input directory	
    os.chdir(input_dir)
   
    #Create component directory 'generate_hcv_consensus' and check if it already exists
    	
    component_dir = "generate_hcv_consensus"
    if not os.path.exists(component_dir):
        os.makedirs(component_dir)
        log_writer.write_log(logger,
             'Created component directory generate_hcv_consensus',
             'info')
    else:
        log_writer.write_log(logger,
             'ERROR: Directory %s exists!' % (component_dir),
             'error')
        sys.exit(67)
     
    #Check for presence of gzipped processed fastq files
    fastqs_exist = len(glob.glob(os.path.join(input_dir, '*.processed.R*.fastq.gz'))) == 2

    if not fastqs_exist:
        log_writer.write_log(logger,
            "ERROR: No processed fastqs in %s" % (input_dir),
            "error")
	sys.exit(68)

    for file in glob.glob(r'*processed*fastq*'):
        print file

        # Copy processed fastq files to generate_hcv_consensus directory for VICUNA
        shutil.copy2(file, component_dir)
        sample_detail = file.split('.')
        sample = sample_detail[0] + '.' + sample_detail[1] + '.' + sample_detail[2]
        res_sample = sample_detail[0]
    print sample
    
    #Change to component directory   	
    os.chdir(component_dir)
    # print 'Current dir', os.getcwd()
    log_writer.write_log(logger,
                         'Current dir: %s\n' % (os.getcwd()),
                         "info")

    #Unzip processed fastqs as required for read counting
    for pfq in glob.glob(r'*processed*fastq*'):
        try:
            p1 = subprocess.Popen(['gunzip', pfq],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE)

            (stdoutdata, stderrdata) = p1.communicate()
            
            if p1.returncode == 0:
                log_writer.write_log(logger,
                         'Files were successfully unzipped: %s\n' % (stdoutdata),
                         "info")
            else:
                log_writer.write_log(logger,
                         'ERROR: Failed to unzip files. Exit code:%s\n. Error message: %s\n' % (p1.returncode, stderrdata),
                         "error")
                sys.exit(69)

        except OSError as e:
            sys.exit("failed to execute program '%s': %s" % (p1, str(e)))

#    output_dir = Args.resxml
    
    #Check assembly directory doesn't already exists and creates one if not
    assembly_dir = 'assembly/'
    if not os.path.exists(assembly_dir):
        os.makedirs(assembly_dir)
        log_writer.write_log(logger,
                         'Assembly dir: %s\n' % (assembly_dir),
                         "info")


#    NEED TO PARSE THIS FILENAME APPROPRIATELY
#      NGS4_Sample_191-1.HCV-sequence-capture.ngsservice.processed.R2.fastq.gz

    #Set paths to various files and refsets    
    hq_fq1 = '%s.processed.R1.fastq' % (sample)
    hq_fq2 = '%s.processed.R2.fastq' % (sample)    
    filt_fq1 = 'assembly/%s_filtered_1.fastq' % (sample)
    filt_fq2 = 'assembly/%s_filtered_2.fastq' % (sample)
    db = '{}/hg38_hcv_nless_k15_s3'.format(Args.refset_dir)
    hcvfasta = '{}/lastz_hcv.fasta'.format(Args.refset_dir)
    contigs= 'assembly/%s.contigs.fasta' % (sample)
    best_ref_fasta = '%s-ref.fasta' % (sample)
    lastz_path = 'lastz'
    lastz_analysed_file = 'assembly/lastz_analysed_file'
    cons1_sorted_final = 'assembly/%s.consensus1.sorted.bam' % (sample)    
    
    human_filtering (sample, db, hq_fq1, hq_fq2, logger, filt_fq1, filt_fq2)
    print 'Filtering', filt_fq1
    split_pops(logger, Args, sample, filt_fq1, filt_fq2)
    print 'Splitpops', filt_fq1
    denovo_assembly(logger, sample, Args)
    print 'De novo', filt_fq1
    find_best_ref(sample, hcvfasta, logger)
    print 'Find bestref', filt_fq1
    assemble_draft(sample, lastz_path, best_ref_fasta, logger)
    print 'Assemble draft', filt_fq1
    contig_map(sample, contigs, filt_fq1, filt_fq2, best_ref_fasta, lastz_analysed_file, logger)
    print 'Contig map', filt_fq1
    quasibam(sample, cons1_sorted_final, logger)
    print 'Quasibam', filt_fq1
    cat_xmls(input_dir, component_dir, logger, res_sample)
    print 'Combine ghc XMLs', sample

    #Write ComponentComplete.txt to signify successful completion of component
    write_component_complete('%s/%s' % (Args.input, 'generate_hcv_consensus'))
   
    return 0

def human_filtering (sample, db, hq_fq1, hq_fq2, logger, filt_fq1, filt_fq2):
    """Requires HCV smalt database index path and trimmed fastq files"""

    '''
    Runs smalt to map reads against a combined human and HCV reference set.
    Human reads are then filtered out for remaining analysis.

    Parameters
    ----------
    sample: str
        sample name
    db: str
        indexed smalt database
    hq_fq1: str
        trimmed forward fastq file
    hq_fq2: str
        trimmed reverse fastq file
    logger: obj
        logger object
    filt_fq1: str
        filtered forward fastq file
    filt_fq2: str
        filtered reverse fastq file

    Returns
    -------
    result: files
        results in 2 human read filtered fastq files for further processing
    '''


    pairs_sam = '%s_pairs.sam' % (sample)

    with open(hq_fq1) as hq1:
        for i, l in enumerate(hq1):
            pass

    fw_reads = (i + 1)/4

    print 'Trimmed forward reads', fw_reads

    with open(hq_fq2) as hq2:
        for i, l in enumerate(hq2):
            pass
    log_writer.write_log(logger,
                        'Forward reads: %s\n' % (fw_reads),
                         "info")
    rev_reads = (i + 1)/4

    print 'Trimmed reverse reads', rev_reads

    log_writer.write_log(logger,
                         'Reverse reads: %s\n' % (rev_reads),
                         "info")

    with open(pairs_sam, 'wb') as sam:

        #Check if sam file exists
        if not os.path.exists(pairs_sam):
            log_writer.write_log(logger,
                             'ERROR: %s file does not exist!\n' % (pairs_sam),
                             "error")
            sys.exit(70)

        if os.path.exists(pairs_sam):

            #Run smalt alignment
        
            p1 = subprocess.Popen(['smalt', 'map', '-x', 
                                   '-y', '0.5', 
                                   '-i', '500', 
                                   '-n', '8', 
                                   db, hq_fq1, hq_fq2],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE)
            
            #Filter reads using awk
            p2 = subprocess.Popen(['awk', """{if ($3 !~ /^chr/ && $7 !~ /^chr/) print $0}"""],            
            #Pass smalt stdout to grep stdin
            stdin=p1.stdout,
            #Pass stdout to sam file
            stdout=sam,
            stderr=subprocess.PIPE)
            #Close file handle for sam file
            p1.stdout.close()
        
            (stdoutdata, stderrdata) = p2.communicate()
    
    for flag, fq in zip(('64','128'), (filt_fq1, filt_fq2)):
        with open(fq, 'wb') as openfq:
            
#            bam = '%s_%s_pairs.bam' % (sample, flag)
      
            #Check if bam file exists
#            if not os.path.exists(bam):
#                log_writer.write_log(logger,
#                         'ERROR: %s file is does not exist!\n' % (bam),
#                         "error")
#                sys.exit(72)

            #Check if pairs sam file produced by smalt is not empty
            if os.stat(pairs_sam).st_size == 0:
                log_writer.write_log(logger,
                         'ERROR: %s file is empty\n' % (pairs_sam),
                         "error")
                sys.exit(71)

#            if os.path.exists(bam):

      
            #Changed phred from forward 64 to 33 and reverse 128 to 33
            #Samtools view to convert SAM to BAM file
            try:
                p1 = subprocess.Popen(['samtools', 'view', '-bhf', flag, pairs_sam],
#               p1 = subprocess.check_call(['samtools', 'view', '-bhf', flag, '-F', '0xC', pairs_sam])
                #Pass stdout to bam file
                #stdout=open(bam,'w'),
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE)
            
                p2 = subprocess.Popen(['samtools', 'bam2fq', '-'],
                stdin=p1.stdout,
                stdout=openfq,
                stderr=subprocess.PIPE)

                (stdoutdata, stderrdata) = p2.communicate()
              
                if p1.returncode == 0:
                    log_writer.write_log(logger,
                                  'Samtools view ran successfully: %s\n' % (stdoutdata),
                                  "info")

		elif '[M::main_bam2fq]' in stderrdata:
                    log_writer.write_log(logger,
                                  'Samtools view ran successfully: %s\n' % (stdoutdata),
                                  "info")
                else:
                    log_writer.write_log(logger,
                                  'ERROR: Samtools view failed. Exit code:%s\n. Error message: %s\n' % (p1.returncode, stderrdata),
                                  "error")
                    sys.exit(72)


                if p2.returncode == 0:
                    log_writer.write_log(logger,
                                  'Samtools bam2fq ran successfully: %s\n' % (stdoutdata),
                                  "info")
                else:
                    log_writer.write_log(logger,
                                  'ERROR: Samtools bam2fq failed. Exit code:%s\n. Error message: %s\n' % (p2.returncode, stderrdata),
                                  "error")
                    sys.exit(73)

            except OSError as e:
                sys.exit("failed to execute program '%s': %s" % (p1, str(e)))
               
    
    # Count the number of reads in the human filtered fastqs

    with open(filt_fq1) as f1:
        for i, l in enumerate(f1):
            pass

    filt_fw_reads = (i + 1)/4

    print 'Filtered forward reads', filt_fw_reads

    log_writer.write_log(logger,
                         'Filtered forward reads: %s\n' % (filt_fw_reads),
                         "info")

    with open(filt_fq2) as f2:
        for i, l in enumerate(f2):
            pass

    filt_rev_reads = (i + 1)/4

    print 'Filtered reverse reads', filt_rev_reads
    
    log_writer.write_log(logger,
                         'Filtered reverse reads: %s\n' % (filt_rev_reads),
                         "info")

    res_type = ET.Element('result')
    res_type.set('type','ghc trimming')
    item1 = ET.SubElement(res_type, 'result_data')
    item2 = ET.SubElement(res_type, 'result_data')
    item3 = ET.SubElement(res_type, 'result_data')
    item4 = ET.SubElement(res_type, 'result_data')
    item1.set('type','fwd trimmed reads')
    item2.set('type','rev trimmed reads')
    item3.set('type','fwd filtered reads')
    item4.set('type','rev filtered reads')
    item1.set('value', str(fw_reads))
    item2.set('value', str(rev_reads))
    item3.set('value', str(filt_fw_reads))
    item4.set('value', str(filt_rev_reads))
    
    xmldata = minidom.parseString(
        ET.tostring(res_type, encoding='utf-8')
        ).toprettyxml(indent="  ")
    samplexml = '%s_trim_filter.xml' % sample
    xmlfile = open(samplexml, 'w')
    xmlfile.write(xmldata)
    
def split_pops(logger, Args, sample, filt_fq1, filt_fq2):
#    ###SPLIT POPULATIONS###
    
    """Processes the filtered fastqs through the snork splitpops software that aims
    to identify presence of mixed infections"""

    '''
    Maps individual reads against a database of HCV references sequences.

    Parameters
    ----------
    sample: str
        sample name
    logger: obj
        logger object
    Args: str
        reference database
    filt_fq1: str
        filtered forward fastq file
    filt_fq2: str
        filtered reverse fastq file

    Returns
    -------
    result: file
        results in a text output of the percentages of each genotype present within the sample
    '''


    #Essential that this directory is completely separate as splitpops cleans up the directory and deletes all input files
    split_dir = 'splitpops/'
    if not os.path.exists(split_dir):
        os.makedirs(split_dir)
        log_writer.write_log(logger,
                         '%s directory created\n' % (split_dir),
                         "info")

    snork = 'snork.py'
    snork_cfg = '{}/snork.config.wtchg'.format(Args.refset_dir)
#    print 'SNORK CONFIG:\n', snork_cfg
    target_ref = '{}/new_hcvrefset'.format(Args.refset_dir)
#    target_ref = '/home/kieren/generate_hcv_consensus/1-dev3/mixed_infection_analysis/new_hcvrefset'   
    filt_bam = '%s_filtered.bam' % (sample)

    #Check if filtered bam file exists
#    if not os.path.exists(bam):
#        log_writer.write_log(logger,
#                         'ERROR: %s file is does not exist!\n' % (filt_bam),
#                         "error")
#        sys.exit(75)

    #Check if filtered bam file is not empty
#    if os.stat(filt_bam).st_size == 0:
#        log_writer.write_log(logger,
#                         'ERROR: %s file is empty\n' % (filt_bam),
#                         "error")
#        sys.exit(76)
#
#    if os.path.exists(filt_bam):


#    #Convert filtered fastqs into bam file
#    filt_bam = '%s_filtered.bam' % (sample)
 #       try:
    p1 = subprocess.Popen(['java', '-jar',
#                           '$PICARD_TOOLS_PATH/FastqToSam.jar',
			   '/phengs/hpc_software/picard-tools/1.111/FastqToSam.jar',	
                           'F1='+ filt_fq1,
                           'F2='+ filt_fq2,
                           'O='+ filt_bam,
                           'SM='+ sample],    
    stdout=subprocess.PIPE, 
    stderr=subprocess.PIPE)


    (stdoutdata, stderrdata) = p1.communicate()
        
    if p1.returncode == 0:
        log_writer.write_log(logger,
                            'FastqToSam ran successfully: %s\n' % (stdoutdata),
                            "info")
    else:
        log_writer.write_log(logger,
                            'ERROR: FastqToSam failed. Exit code:%s\n. Error message: %s\n' % (p1.returncode, stderrdata),
                            "error")
        sys.exit(74)
        
    #Run splitpops.py on bam file

#    p2 = subprocess.list2cmdline([splitpops, 'splitpops', '-bin snorktest/src6',
#        try:

        #Check if filtered bam file exists
    if not os.path.exists(filt_bam):
        log_writer.write_log(logger,
                             'ERROR: %s file is does not exist!\n' % (filt_bam),
                             "error")
        sys.exit(75)

    #Check if filtered bam file is not empty
    if os.stat(filt_bam).st_size == 0:
        log_writer.write_log(logger,
                             'ERROR: %s file is empty\n' % (filt_bam),
                             "error")
        sys.exit(76)

    if os.path.exists(filt_bam):

        p2 = subprocess.Popen([snork, 'splitpops',
                               '-bin', 'snorktest/src6',
                               '-profile', 'None',
                               '-config', snork_cfg,
                               '-orgid', 'Hepc',
                               '-dataid', sample,
                               '-samplename', sample,
                               '-bampath', filt_bam,
                               '-targetrefid', 'new_hcvrefset',
                               '-targetrefpath', target_ref,
                               '-outdir', split_dir,
                               '-logdir', split_dir,
                               '-overwrite', 'False',
                               '-deleteints', 'True',
                               '-verbosity', 'DEBUG'],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE)
        
        (stdoutdata, stderrdata) = p2.communicate()

        if p2.returncode == 0:
            log_writer.write_log(logger,
                                 'Snork splitpops ran successfully: %s\n' % (stdoutdata),
                                 "info")
        else:
            log_writer.write_log(logger,
                                 'ERROR: Snork splitpops failed. Exit code:%s\n. Error message: %s\n' % (p2.returncode, stderrdata),
                                 "error")
            sys.exit(77)


    split_stats = 'splitpops/%s_filtered.targetpop.stats.txt' % (sample)
    
    #Create dictionary of genotypes(k) and percentage reads(v) and return the highest two
    d={}
    with open (split_stats, 'r') as s:
        #Skips the header line and only returns the top 2 lines (Highest scoring genotype matches)
        lines = s.readlines()[1:3]
        main_geno = lines[0].split()
        main_geno_name = main_geno[0]
        main_geno_perc = main_geno[2]
        second_geno = lines[1].split()
        second_geno_name = second_geno[0]
        second_geno_perc = second_geno[2]

	#Create splitpops XML
        res_type = ET.Element('result')
        res_type.set('type','ghc genotypes')
        item1 = ET.SubElement(res_type, 'result_data')
        item2 = ET.SubElement(res_type, 'result_data')
        item1.set('type', main_geno_name)
        item1.set('value', str(main_geno_perc))
        item2.set('type', second_geno_name)
        item2.set('value', str(second_geno_perc))

        xmldata = minidom.parseString(
            ET.tostring(res_type, encoding='utf-8')
            ).toprettyxml(indent="  ")
        samplexml = '%s_splitpops.xml' % sample
        xmlfile = open(samplexml, 'w')
        xmlfile.write(xmldata)


#    #Run bamtools split on splitpops bam
#    bamtools = 'bamtools split -in "temp/"$sample"_filtered.targetpop.readclass.bam" -tag RG:Z ;'
#    os.system(bamtools)
#    #Select largest file from using *RG* tag
#    largest_bam = '`find temp/*RG* -type f | xargs ls -1S | head -n 1`' #Need to check
#    os.system(largest_bam)
#    #Read file for largest bam file filename
#    #largest_bam= cat "temp/"$sample"_splitpop_largest.txt" ;
#
#    #Sort largest bam by readname and then convert selected largest splitpop bam to fastqs
#    sort_largest_bam = 'samtools sort -n $largest_bam_path $largest_bam_path"_sorted"'
#    os.system(sort_largest_bam)
#    bam_to_fastq = 'bamToFastq -i $largest_bam_path"_sorted.bam" -fq "temp/"$sample"_splitpop_1.fastq" -fq2 "temp/"$sample"_splitpop_2.fastq"'
#

def denovo_assembly(logger, sample, Args):
    """De novo assembly using VICUNA accesses a preconfigured config file that
    contains the paths to the directory containing the fastq files."""

    '''

    Parameters
    ----------
    sample: str
        sample name
    logger: obj
        logger object
    Args: str
        reference set
    Returns
    -------
    result: file
        results in a file containing all the generated contigs
    '''

    ###Starting VICUNA###
    #Run VICUNA on sample specific config file
    try:      
        vicuna_config = os.environ['GENERATE_HCV_CONSENSUS_VICUNA_CONFIG']#'/phengs/hpc_software/ucl_assembly/vicuna_config.txt'
    except KeyError, e:
        log_writer.write_log(logger,
                            'Error: GENERATE_HCV_CONSENSUS_VICUNA_CONFIG environment variable not set: %s\n' % (sys.exc_info()[0]),
                            "error")
        sys.exit(78)
#        raise EnvironmentError('GENERATE_HCV_CONSENSUS_VICUNA_CONFIG environment variable needs to be set to import this module') from KeyError
    
    with open(vicuna_config) as vc:
        cwd = os.getcwd()
        s = vc.read().replace('*DIR*', cwd).replace('*SAMPLE*', sample).replace('*MSA_FILE*', '{}/hcv_wgs_alignment.fas'.format(Args.refset_dir))

    vicuna_config_sample = 'vicuna_config_%s.txt' % (sample)
    with open(vicuna_config_sample, 'w') as vc:
        vc.write(s)

    #Create directory for VICUNA assembly and check to see if it already exists
#    vicuna_dir_path = '%s/assembly/' % (sample)
#    vicuna_dir = os.path.dirname(vicuna_dir_path)
#    if not os.path.exists(vicuna_dir):
#               os.makedirs(vicuna_dir)
    
    #Run VICUNA with the new sample config file

    p1 = subprocess.Popen(['vicuna', vicuna_config_sample],
   
    stdout=subprocess.PIPE,
    stderr=subprocess.PIPE)

    (stdoutdata, stderrdata) = p1.communicate()

    if p1.returncode == 0:
        log_writer.write_log(logger,
                             'VICUNA ran successfully: %s\n' % (stdoutdata),
                             "info")
    else:
        log_writer.write_log(logger,
                             'ERROR: VICUNA failed. Exit code:%s\n. Error message: %s\n' % (p1.returncode, stderrdata),
                             "error")
        sys.exit(79)

    #Following completion of VICUNA, the 'dg-0' needs to be replaced by sample name and 
    # contig number starting at 1 as VICUNA starts at 0
    # and each contig needs a specific whole number for the hash later in the pipeline

    contigs = 'assembly/contig.fasta'
    newdata = []    

    #Open contigs file for reading and replacing header line with new information on sample and sensible contig number
    with open(contigs) as c:
        for line in c.readlines():
            if line.startswith('>dg-'):
                start_line = re.search("^\S*", line).group(0)
                element = start_line.split('-')

                #Returns just the number from the match object reflecting the VICUNA contig number
                num = element[1]
                print num

                #Adds an increment of 1 to indicate new contig numbering
                count = int(num) + 1
                #Number of contigs
                count = str(count)
                newdata.append(re.sub(r'dg-.*', sample + '.' + count, line).strip())
                #print newdata
            else:
                #print line.strip()        
                newdata.append(line.strip())
    
    #Create denovo XML
    res_type = ET.Element('result')

    if int(count) == 0:
        print 'No contig'
        res_type.set('type','ghc denovo')
        item1 = ET.SubElement(res_type, 'result_data')
        item1.set('type','Contig number')
        item1.set('value','No contigs')


    elif int(count) == 1:
        print 'Single contig'
        res_type.set('type','ghc denovo')
        item2 = ET.SubElement(res_type, 'result_data')
        item2.set('type','Contig number')
        item2.set('value','Single contig')


    else:
        print 'Number of contigs =', count
        res_type.set('type','ghc denovo')
        item3 = ET.SubElement(res_type, 'result_data')
        item3.set('type','Contig number')
        item3.set('value', count)


    print '\n'.join(newdata)

    
    # Print pretty XML  
    xmldata = minidom.parseString(
        ET.tostring(res_type, encoding='utf-8')
        ).toprettyxml(indent="  ")   
#    xmldata = ET.tostring(result)
    samplexml = '%s_denovo.xml' % sample
    xmlfile = open(samplexml, 'w')
    xmlfile.write(xmldata)     
    
    #Very important that this file is named correctly as lastz_analyser.WITH_REVCOMP.pl script will fail later in the pipeline.
    # It must be named samplename.contigs.fasta
    contigs_increment = 'assembly/%s.contigs.fasta' % (sample)
    
    with open(contigs_increment, 'w') as ci:
        ci.write('\n'.join(newdata))

    #Remove processed fastq copies in the assembly directory
    for file in glob.glob(r'*processed*fastq*'):
        os.remove(file)
                
def find_best_ref(sample, hcvfasta, logger):
    """Following completion of de novo assembly, the best reference sequence from the HCV database needs to
    be selected for mapping to create a draft assembly"""

    '''
    Runs lastz to map contigs against the HCV reference set.

    Parameters
    ----------
    sample: str
        sample name
    hcvfasta: str
        HCV reference set
    logger: obj
        logger object

    Returns
    -------
    result: files
        returns the optimum HCV reference sequence
    '''

    
    #Path to current available lastz executable
    lastz_path = 'lastz'
    #Path to contigs fasta file with the specification of [multiple] required by lastz if multiple contigs
    contigs = 'assembly/%s.contigs.fasta[multiple]' % (sample)
    #Path to lastz output file
    contigs_lastz = 'assembly/contig.lastz'
    
    
    with open(contigs_lastz, 'w') as cz:
            #Run LASTZ
    
        p1 = subprocess.Popen([lastz_path, contigs, hcvfasta, 
                              '--ambiguous=iupac',
                              '--format=GENERAL'],
        #Pass stdout to contigs lastz file
        stdout=cz,
        stderr=subprocess.PIPE)
    #   p1.stdout.close()

        (stdoutdata, stderrdata) = p1.communicate()

        if p1.returncode == 0:
            log_writer.write_log(logger,
                                 'Lastz ran successfully: %s\n' % (stdoutdata),
                                 "info")
        else:
            log_writer.write_log(logger,
                                 'ERROR: Lastz failed. Exit code:%s\n. Error message: %s\n' % (p1.returncode, stderrdata),
                                 "error")
            sys.exit(80)

    lastz_log = 'assembly/lastz_besthit.log'
    lastz_bestref = 'lastz_bestref.pl'
    best_ref_fasta = '%s-ref.fasta' % (sample)

    #Run lastz_bestref.pl to find the best matching reference sequence. The perl script requires string concatenation between the flags
    # and arguments (+) for these to run

    p1 = subprocess.Popen([lastz_bestref, 
                           '-contig_lastz='+ contigs_lastz,
                           '-blastdb='+ hcvfasta,
                           '-best_ref_fasta='+ best_ref_fasta,
                           '-lastz_best_hit_log='+ lastz_log],
    stdout=subprocess.PIPE,
    stderr=subprocess.PIPE)
        
    (stdoutdata, stderrdata) = p1.communicate()
    
    if p1.returncode == 0:
        log_writer.write_log(logger,
                             'Lastz bestref ran successfully: %s\n' % (stdoutdata),
                             "info")
    else:
        log_writer.write_log(logger,
                             'ERROR: Lastz bestref failed. Exit code:%s\n. Error message: %s\n' % (p1.returncode, stderrdata),
                             "error")
        sys.exit(81)

def assemble_draft(sample, lastz_path, best_ref_fasta, logger):
    """The best reference has been defined, now assemble the draft genome"""

    '''

    Parameters
    ----------
    sample: str
        sample name
    lastz_path: str
        lastz path
    logger: obj
        logger object
    best_ref_fasta: str
        Optimum reference sequence

    Returns
    -------
    result: files
        Returns assembled draft sequence
    '''


    contigs = 'assembly/%s.contigs.fasta' % (sample)
    contigs_bestref = 'assembly/contigs-vs-bestref.lav'


    with open(contigs_bestref, 'w') as cb:
    
        #Run lastz to compare contigs to best reference. NOTE: NEED TO ENSURE THE FUNCTION ARGUMENTS ARE PASSED FROM PREVIOUS CORRECTLY
        p1 = subprocess.Popen([lastz_path, best_ref_fasta, contigs,
                               '--ambiguous=iupac'],
        #Pass stdout to contigs_vs_bestref.lav file
        stdout=cb,
        stderr=subprocess.PIPE)

        (stdoutdata, stderrdata) = p1.communicate()
        
        if p1.returncode == 0:
            log_writer.write_log(logger,
                                 'Lastz compare ran successfully: %s\n' % (stdoutdata),
                                 "info")
        else:
            log_writer.write_log(logger,
                                 'ERROR: Lastz compare failed. Exit code:%s\n. Error message: %s\n' % (p1.returncode, stderrdata),
                                 "error")
            sys.exit(82)

    lastz_analyser_rev = 'lastz_analyser.WITH_REVCOMP.pl'
    lastz_analysed_file = 'assembly/lastz_analysed_file'
    lastz_analysed_log = 'assembly/lastz_analyser.log'

 
    p2 = subprocess.Popen([lastz_analyser_rev,
                           '-reference_fasta_file='+ best_ref_fasta,
                           '-sample_fasta_file='+ contigs,
                           '-lastz_results_file='+ contigs_bestref,
                           '-cutoff=50000',
                           '-with_revcomp=yes',
                           '-output='+ lastz_analysed_file,
                           '-log_file='+ lastz_analysed_log],
    #Pass stdout to contigs lastz file
    stdout=subprocess.PIPE,
    stderr=subprocess.PIPE)

    (stdoutdata, stderrdata) = p2.communicate()
    
    if p2.returncode == 0:
        log_writer.write_log(logger,
                             'Lastz analyser reverse ran successfully: %s\n' % (stdoutdata),
                             "info")
    else:
        log_writer.write_log(logger,
                             'ERROR: Lastz analyser reverse failed. Exit code:%s\n. Error message: %s\n' % (p2.returncode, stderrdata),
                             "error")
        sys.exit(83)

def contig_map(sample, contigs, filt_fq1, filt_fq2, best_ref_fasta, lastz_analysed_file, logger):
    """The best reference has been defined, now assemble the draft genome"""

    '''

    Parameters
    ----------
    sample: str
        sample name
    contigs: str
        contigs file
    filt_fq1: str
        filtered forward fastq file
    filt_fq2: str
        filtered reverse fastq file
    best_ref_fasta: str
        Optimum reference sequence
    lastz_analysed_file: str
        lastz result file
    logger: obj
        logger object

    Returns
    -------
    result: file
        Returns consensus sequence
    '''

    ## map reads using BWA MEM to contigs to work out which sequences to choose when contigs overlap

    contigs_sam = 'assembly/%s.contigs.sam' % (sample)
    contigs_bam = 'assembly/%s.contigs.bam' % (sample)
    contigs_sorted = 'assembly/%s.contigs.sorted' % (sample)
    contigs_sorted_final = 'assembly/%s.contigs.sorted.bam' % (sample)
    contigs_mpileup = 'assembly/%s.contigs.mpileup' % (sample)
    genome_maker = 'genome_maker2b.pl'
    genome_fasta = 'assembly/%s-genome.fasta' % (sample)
    genome_maker_log = 'assembly/%s_genome_maker.log' % (sample)
    genome_sam = 'assembly/%s-genome.sam' % (sample)
    genome_bam = 'assembly/%s-genome.bam' % (sample)
    genome_sorted = 'assembly/%s-genome.sorted' % (sample)
    genome_sorted_final = 'assembly/%s-genome.sorted.bam' % (sample)
    genome_mpileup = 'assembly/%s-genome.mpileup' % (sample)
    cons_mv = 'cons_mv.pl'
    cons_pren = 'assembly/%s.consensus1.preNcut.fasta' % (sample)
    genome_mv = 'assembly/%s-genome.fasta.mv' % (sample)
    genome_basefreq = 'assembly/%s-genome.fasta.basefreqs.tsv' % (sample)
    n_remover = 'N_remover_from_consensus.pl'
    cons1_fasta =    'assembly/%s.consensus1.fasta' % (sample)
    cons1_sam = 'assembly/%s.consensus1.sam' % (sample)
    cons1_bam = 'assembly/%s.consensus1.bam' % (sample)
    cons1_sorted = 'assembly/%s.consensus1.sorted' % (sample)
    cons1_sorted_final = 'assembly/%s.consensus1.sorted.bam' % (sample)
    cons1_mpileup = 'assembly/%s.consensus1.mpileup' % (sample)
    cons2_pren = 'assembly/%s.consensus2.preNcut.fasta' % (sample)
    cons1_mv = 'assembly/%s.consensus1.fasta.mv' % (sample)
    cons1_basefreq = 'assembly/%s.consensus1.fasta.basefreqs.tsv' % (sample)
    cons2_fasta = 'assembly/%s.consensus2.fasta' % (sample)
    majvar = 'majvarcheck2_bwa.pl'

    cmd_gm = [genome_maker,
              '-sample_pileup_file='+ contigs_mpileup,
              '-contigs='+ contigs,
              '-reference_mapped_consensus='+ best_ref_fasta,
              '-lastz_analysed_file='+ lastz_analysed_file,
              '-ref_correct_start=0',
              '-ref_correct_stop=20000',
              '-output='+ genome_fasta,
              '-logfile='+ genome_maker_log]

    for fasta, sam, bam, sorted, final, mpileup in zip((contigs, genome_fasta, cons1_fasta),
                                                       (contigs_sam, genome_sam, cons1_sam),
                                                       (contigs_bam, genome_bam, cons1_bam),
                                                       (contigs_sorted, genome_sorted, cons1_sorted),
                                                       (contigs_sorted_final, genome_sorted_final, cons1_sorted_final),
                                                       (contigs_mpileup, genome_mpileup, cons1_mpileup)):
    
        p1 = subprocess.Popen(['bwa', 'index', fasta],
        
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE)
    
        (stdoutdata, stderrdata) = p1.communicate()

        if p1.returncode == 0:
            log_writer.write_log(logger,
                                'BWA index ran successfully: %s\n' % (stdoutdata),
                                "info")
        else:
            log_writer.write_log(logger,
                                 'ERROR: BWA index failed. Exit code:%s\n. Error message: %s\n' % (p2.returncode, stderrdata),
                                 "error")
            sys.exit(84)


        with open(sam, 'w') as s:

        
            p1 = subprocess.Popen(['bwa', 'mem',
                                   '-t', '8', 
                                   fasta, 
                                   filt_fq1, 
                                   filt_fq2],
            stdout=s,
            stderr=subprocess.PIPE)
    
            (stdoutdata, stderrdata) = p1.communicate()
            
            if p1.returncode == 0:
                log_writer.write_log(logger,
                                     'BWA MEM ran successfully: %s\n' % (stdoutdata),
                                     "info")
            else:
                log_writer.write_log(logger,
                                     'ERROR: BWA MEM failed. Exit code:%s\n. Error message: %s\n' % (p1.returncode, stderrdata),
                                     "error")
                sys.exit(85)
            

            p2 = subprocess.Popen(['samtools', 'view', '-bhS',
                                   '-o', bam, sam],
            stdin=p1.stdout,    
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE)
    
            (stdoutdata, stderrdata) = p2.communicate()
            
            if p2.returncode == 0:
                log_writer.write_log(logger,
                                     'Samtools view ran successfully: %s\n' % (stdoutdata),
                                     "info")
            else:
                log_writer.write_log(logger,
                                     'ERROR: Samtools view failed. Exit code:%s\n. Error message: %s\n' % (p2.returncode, stderrdata),
                                     "error")
                sys.exit(86)
    
            p3 = subprocess.Popen(['samtools', 'sort',
                                   '-@', '8', bam, sorted],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE)

            (stdoutdata, stderrdata) = p3.communicate()

            if p3.returncode == 0:
                log_writer.write_log(logger,
                                     'Samtools sort ran successfully: %s\n' % (stdoutdata),
                                     "info")
            else:
                log_writer.write_log(logger,
                                    'ERROR: Samtools sort failed. Exit code:%s\n. Error message: %s\n' % (p3.returncode, stderrdata),
                                    "error")
                sys.exit(87)

        
            p4 = subprocess.Popen(['samtools', 'index', final],
   
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE)

            (stdoutdata, stderrdata) = p4.communicate()

            if p4.returncode == 0:
                log_writer.write_log(logger,
                                     'Samtools index ran successfully: %s\n' % (stdoutdata),
                                     "info")
            else:
                log_writer.write_log(logger,
                                    'ERROR: Samtools index failed. Exit code:%s\n. Error message: %s\n' % (p4.returncode, stderrdata),
                                    "error")
                sys.exit(88)
       
 
            p5 = subprocess.Popen(['samtools', 'mpileup', '-f',
                                   fasta,
                                   '-d', '1000000',
                                   '-o', mpileup, final],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE)

            (stdoutdata, stderrdata) = p5.communicate()

            if p5.returncode == 0:
                log_writer.write_log(logger,
                                     'Samtools mpileup ran successfully: %s\n' % (stdoutdata),
                                     "info")
            else:
                log_writer.write_log(logger,
                                    'ERROR: Samtools mpileup failed. Exit code:%s\n. Error message: %s\n' % (p5.returncode, stderrdata),
                                    "error")
                sys.exit(89)

    
        #Helpful method for checking exactly what is sent to the command line (subprocess.list2cmdline)
            if fasta==contigs:    

#               p6 = subprocess.list2cmdline(cmd_gm)
                p6 = subprocess.Popen(cmd_gm, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    
                #Print the full command line argument
                print 'fasta is contigs'

                (stdoutdata, stderrdata) = p6.communicate()

                if p6.returncode == 0:
                    log_writer.write_log(logger,
                                         'Genome maker ran successfully: %s\n' % (stdoutdata),
                                         "info")
                else:
                    log_writer.write_log(logger,
                                        'ERROR: Genome maker failed. Exit code:%s\n. Error message: %s\n' % (p6.returncode, stderrdata),
                                        "error")
                    sys.exit(90)


            elif fasta==genome_fasta:

                p7 = subprocess.Popen([cons_mv,
                                      '-mpileup='+ genome_mpileup,
                                      '-reference_fasta='+ genome_fasta,
                                      '-mv_freq_cutoff=0.01',
                                      '-mv_variant_depth_cutoff=20',
                                      '-cons_depth_cutoff=80',
                                      '-sliding_window_size=300',
                                      '-consensus_out='+ cons_pren,
                                      '-mv_out='+ genome_mv,
                                      '-base_freq_out='+ genome_basefreq], 
                
                stdout=subprocess.PIPE, 
                stderr=subprocess.PIPE)

                print 'fasta is genome'
                (stdoutdata, stderrdata) = p7.communicate()

                if p7.returncode == 0:
                    log_writer.write_log(logger,
                                         'Cons mv ran successfully: %s\n' % (stdoutdata),
                                         "info")
                else:
                    log_writer.write_log(logger,
                                         'ERROR: Cons mv failed. Exit code:%s\n. Error message: %s\n' % (p7.returncode, stderrdata),
                                         "error")
                    sys.exit(91)
               
 
                p8 = subprocess.Popen([n_remover,
                                      '-cutoff=46',
                                      cons_pren],

                stdout=open(cons1_fasta, 'w'), 
                stderr=subprocess.PIPE)                
                
                (stdoutdata, stderrdata) = p8.communicate()

                if p8.returncode == 0:
                    log_writer.write_log(logger,
                                         'N Remover ran successfully: %s\n' % (stdoutdata),
                                         "info")
                else:
                    log_writer.write_log(logger,
                                        'ERROR: N Remover failed. Exit code:%s\n. Error message: %s\n' % (p8.returncode, stderrdata),
                                         "error")
                    sys.exit(92)


            else:

                p9 = subprocess.Popen([cons_mv,
                                      '-mpileup='+ cons1_mpileup,
                                      '-reference_fasta='+ cons1_fasta,
                                      '-mv_freq_cutoff=0.01',
                                      '-mv_overall_depth_cutoff=100',
                                      '-mv_variant_depth_cutoff=20',
                                      '-cons_depth_cutoff=80',
                                      '-sliding_window_size=300',
                                      '-consensus_out='+ cons2_pren,
                                      '-mv_out='+ cons1_mv,
                                      '-base_freq_out='+ cons1_basefreq],

                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE)

                print 'fasta is consensus'
                (stdoutdata, stderrdata) = p9.communicate()

                if p9.returncode == 0:
                    log_writer.write_log(logger,
                                         'Cons mv basefreq ran successfully: %s\n' % (stdoutdata),
                                         "info")
                else:
                    log_writer.write_log(logger,
                                         'ERROR: Cons mv basefreq failed. Exit code:%s\n. Error message: %s\n' % (p9.returncode, stderrdata),
                                         "error")
                    sys.exit(93)


                p10 = subprocess.Popen([n_remover,
                                       '-cutoff=46',
                                        cons2_pren],

                stdout=open(cons2_fasta, 'w'),
                stderr=subprocess.PIPE)

                (stdoutdata, stderrdata) = p10.communicate()

                if p10.returncode == 0:
                    log_writer.write_log(logger,
                                         'N Remover cons2 ran successfully: %s\n' % (stdoutdata),
                                         "info")
                else:
                    log_writer.write_log(logger,
                                         'ERROR: N Remover cons2 failed. Exit code:%s\n. Error message: %s\n' % (p10.returncode, stderrdata),
                                         "error")
                    sys.exit(94)
    
   
                p11 = subprocess.Popen([majvar,
                                       '-mvpath='+ cons1_mv,
                                       '-basefreq='+ cons1_basefreq,
                                       '-fwdreads='+ filt_fq1,
                                       '-revreads='+ filt_fq2], 
                stdout=subprocess.PIPE, 
                stderr=subprocess.PIPE)
                
                (stdoutdata, stderrdata) = p11.communicate()

                if p11.returncode == 0:
                    log_writer.write_log(logger,
                                         'Majvar bwa ran successfully: %s\n' % (stdoutdata),
                                         "info")
                else:
                    log_writer.write_log(logger,
                                         'ERROR: Majvar bwa failed. Exit code:%s\n. Error message: %s\n' % (p11.returncode, stderrdata),
                                         "error")
                    sys.exit(95)
               
 
                p12 = subprocess.Popen(['samtools', 'view',
                                        '-S',
                                        '-F0x4',
                                        '-c', contigs_sam],
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE)

                (stdoutdata, stderrdata) = p12.communicate()

                if p12.returncode == 0:
                    log_writer.write_log(logger,
                                         'Samtools view final ran successfully: %s\n' % (stdoutdata),
                                         "info")
                else:
                    log_writer.write_log(logger,
                                         'ERROR: Samtools view final failed. Exit code:%s\n. Error message: %s\n' % (p12.returncode, stderrdata),
                                         "error")
                    sys.exit(96)

#   Generate mapping XML
    res_type = ET.Element('result')
    res_type.set('type','ghc mapping')
    item1 = ET.SubElement(res_type, 'result_data')
    item1.set('type','reads mapped to contigs')
    item1.set('value', str(stdoutdata).rstrip())

    xmldata = minidom.parseString(
        ET.tostring(res_type, encoding='utf-8')
        ).toprettyxml(indent="  ")
    samplexml = '%s_mapping.xml' % sample
    xmlfile = open(samplexml, 'w')
    xmlfile.write(xmldata)

def quasibam(sample, cons1_sorted_final, logger):
    """Quasibam requires the consensus sequence and consensus bam file to
    produce a minority variants report"""

    '''

    Parameters
    ----------
    sample: str
        sample name
    cons1_sorted_final: str
        consensus bam file
    logger: obj
        logger object

    Returns
    -------
    result: files
        Returns consensus sequence and minority variants report
    '''

#    """Run quasibam on consensus bam file and consensus fasta sequence"""

    consensus = 'assembly/%s.consensus2.fasta' % (sample)    
#    quasifas = '%s/assembly/%s_vicuna_bwa_quasibam.fas' % (sample, sample)
#    quasitxt = '%s/assembly/%s_vicuna_bwa_quasibam.txt' % (sample, sample)
#    quasierr = '%s/assembly/%s_vicuna_bwa_quasibam.err' % (sample, sample)

#    with open(quasifas, 'w') as qf:

#    p1 = subprocess.list2cmdline(['quasi_bam',
#                                       cons1_sorted_final,
#                                       consensus])
#    print p1
    
    p1 = subprocess.Popen(['quasi_bam',
                           cons1_sorted_final, 
                           consensus,
                           '-c', '15'],
#                              '-o1', qf],
#                              '-o2', quasitxt,
#                              '-o3', quasierr],
                
    stdout=subprocess.PIPE,
    stderr=subprocess.PIPE)

    (stdoutdata, stderrdata) = p1.communicate()

    if p1.returncode == 0:
        log_writer.write_log(logger,
                             'Quasibam ran successfully: %s\n' % (stdoutdata),
                             "info")
    else:
        log_writer.write_log(logger,
                             'ERROR: Quasibam failed. Exit code:%s\n. Error message: %s\n' % (p1.returncode, stderrdata),
                             "error")
        sys.exit(97)

def cat_xmls(input_dir, component_dir, logger, res_sample):
    """Concatenates all the XMLs into one document"""

    '''

    Parameters
    ----------
    input_dir: str
        input directory containing all XML files
    component_dir: str
        component directory
    res_sample: str
        sample name
    logger: obj
        logger object

    Returns
    -------
    result: file
        Returns combined XML file
    '''


    xmls_exist = len(glob.glob(os.path.join(input_dir + '/' + component_dir + '/' + '*.xml'))) == 4

    if not xmls_exist:
        log_writer.write_log(logger,
                             "ERROR: Expecting 4 xml files in %s" % (input_dir),
                             "error")
        sys.exit(98)

    with open(input_dir + '/' + component_dir + '/' + res_sample + '.results.xml', 'w') as outfile:
        outfile.write('<ngs_sample id="' +  res_sample + '">\n <workflow value="generate_hcv_consensus" version="ngsservice"/>\n')
        for file in glob.glob(os.path.join('*.xml')):
            with open(file) as f:
                for line in f:
                    if "<?xml" not in line:
                        outfile.write(line)

        outfile.write('</ngs_sample>')


if __name__ == '__main__':

    sys.exit(main())
