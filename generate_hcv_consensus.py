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
import log_writer
import csv
import shutil
from xml.dom import minidom
import xml.etree.ElementTree as ET

__version__ = '0.1'
__date__ = '07Jun2018'
__author__ = 'kieren.lythgow@phe.gov.uk'

def parse_args():
	"""
#	Parge arguments
#	Parameters
#	----------
#	no inputs
#	Returns
#	-------
#	Args: obj
#	arguments object
#	"""

	Description = 'version %s, date %s, author %s' %(__version__, __date__, __author__)

	Parser = argparse.ArgumentParser(description=Description)

	Parser.add_argument('-i', '--input',
                          metavar='STRING',
                          dest='input',
                          required=True,
                          help='Input folder incl. path. [Required.]')

# 	Parser.add_argument('-r', '--resultxml',
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

	Args = parse_args()

	#set up logger
	logger = log_writer.setup_logger('logs/contamination_check_viral.stdout',
                                     'logs/contamination_check_viral.stderr')
    log_writer.info_header(logger, "starting contamination check")
    log_writer.write_log(logger,
                         "processing samples in %s ...\n" % oArgs.input,
                         "info")




#	sample = '171009_8' #define sample
#	sample = sys.argv[1]
#	sampl = Args.input
	input_dir = Args.input
#	sample_dir = os.path.dirname(input_dir)
#	sample_detail = os.path.basename(sample_dir).rsplit('_', 4)
#	sample = sample_detail[0]

#	sample = dir_detail[1] + '_' + dir_detail[2]

	os.chdir(input_dir)

	component_dir = "generate_hcv_consensus"
	os.makedirs(component_dir)

	for file in glob.glob(r'*processed*fastq*'):
    		print file
		shutil.copy2(file, component_dir)
		sample_detail = file.split('.')
		sample = sample_detail[0] + '.' + sample_detail[1] + '.' + sample_detail[2]
		
	print sample
	os.chdir(component_dir)
        print 'Current dir', os.getcwd()

	for pfq in glob.glob(r'*processed*fastq*'):
		p1 = subprocess.Popen(['gunzip', file],
		stdout=subprocess.PIPE,
		stderr=subprocess.PIPE)

#       print p1
		(stdoutdata, stderrdata) = p1.communicate()


#	output_dir = Args.resxml


        assembly_dir = 'assembly/'
#       if not os.path.exists(split_dir):
        os.makedirs(assembly_dir)
        print 'Assembly dir', assembly_dir

#	hq_fq1 = '%s_hq_1.fastq' % (sample)
#	hq_fq2 = '%s_hq_2.fastq' % (sample)

#	NEED TO PARSE THIS FILENAME APPROPRIATELY
#  	NGS4_Sample_191-1.HCV-sequence-capture.ngsservice.processed.R2.fastq.gz
	
	hq_fq1 = '%s.processed.R1.fastq' % (sample)
	hq_fq2 = '%s.processed.R2.fastq' % (sample)	
	filt_fq1 = 'assembly/%s_filtered_1.fastq' % (sample)
	filt_fq2 = 'assembly/%s_filtered_2.fastq' % (sample)
	db = '/home/kieren/result_dir/hg38_hcv_k15_s3' #define HCV db filepath
	hcvfasta = '/home/kieren/UCL_IVA/Data/hcv.fasta'
	contigs= 'assembly/%s.contigs.fasta' % (sample)
	best_ref_fasta = '%s-ref.fasta' % (sample)
	lastz_path = '/home/kieren/lastz-distrib/bin/lastz'
	lastz_analysed_file = 'assembly/lastz_analysed_file'
	cons1_sorted_final = 'assembly/%s.consensus1.sorted.bam' % (sample)	
	
	human_filtering (sample, db, hq_fq1, hq_fq2, filt_fq1, filt_fq2)
	print 'Filtering', filt_fq1
#	split_pops(sample, filt_fq1, filt_fq2)
#	print 'Splitpops', filt_fq1
#	denovo_assembly(sample)
#	print 'De novo', filt_fq1
#	find_best_ref(sample, hcvfasta)
#	print 'Find bestref', filt_fq1
#	assemble_draft(sample, lastz_path, best_ref_fasta)
#	print 'Assemble draft', filt_fq1
#	contig_map(sample, contigs, filt_fq1, filt_fq2, best_ref_fasta, lastz_analysed_file)
#	print 'Contig map', filt_fq1
#	quasibam(sample, cons1_sorted_final)
#	print 'Quasibam', filt_fq1

def human_filtering (sample, db, hq_fq1, hq_fq2, filt_fq1, filt_fq2):
	"""Requires HCV smalt database index path and trimmed fastq files"""

	pairs_sam = '%s_pairs.sam' % (sample)
#	filt_fq1 = 'assembly/%s_filtered_1.fastq' % (sample)
#	filt_fq2 = 'assembly/%s_filtered_2.fastq' % (sample)	i
	
	with open(hq_fq1) as hq1:
		for i, l in enumerate(hq1):
			pass

	fw_reads = (i + 1)/4

	print 'Trimmed forward reads', fw_reads

	with open(hq_fq2) as hq2:
                for i, l in enumerate(hq2):
                        pass

        rev_reads = (i + 1)/4

	print 'Trimmed reverse reads', rev_reads

	with open(pairs_sam, 'wb') as sam:
    		#Run smalt alignment
		try:
			p1 = subprocess.Popen(['smalt', 'map', '-x', 
						'-y', '0.5', 
						'-i', '500', 
						'-n', '8', 
						db, hq_fq1, hq_fq2],
    			stdout=subprocess.PIPE,
    			stderr=subprocess.PIPE)
    		
			#Filter reads using awk
			p2 = subprocess.Popen(['awk', """{if ($3 !~ /^chr/ && $7 !~ /^chr/) print $0}"""],    		
			#Filter reads using grep	
#    			p2 = subprocess.Popen(['grep', '-v', 'SN:chr'],
    			#Pass smalt stdout to grep stdin
			stdin=p1.stdout,
			#Pass stdout to sam file
    			stdout=sam,
    			stderr=subprocess.PIPE)
			#Close file handle for sam file
			p1.stdout.close()
    	
	   	 	(stdoutdata, stderrdata) = p2.communicate()
			if stdoutdata:
				print 'smaltrc ', p1.returncode
				print 'greprc ', p2.returncode
				print "SUCCESS ", stdoutdata
			if stderrdata:
				print 'smaltrc ', p1.returncode
                        	print 'greprc ', p2.returncode
				print "ERROR ",stderrdata.strip()
	
		except OSError as e:
    			print "OSError > ",e.errno
    			print "OSError > ",e.strerror
    			print "OSError > ",e.filename
		except:
    			print "Error > ",sys.exc_info()[0]	
	

	for flag, fq in zip(('64','128'), (filt_fq1, filt_fq2)):
		with open(fq, 'wb') as openfq:
			bam = '%s_%s_pairs.bam' % (sample, flag)
    			#Changed phred from forward 64 to 33 and reverse 128 to 33
			#Samtools view to convert SAM to BAM file
			try:
				p1 = subprocess.Popen(['samtools', 'view', '-bhf', flag, pairs_sam],
#				p1 = subprocess.check_call(['samtools', 'view', '-bhf', flag, '-F', '0xC', pairs_sam])
				#Pass stdout to bam file
				#stdout=open(bam,'w'),
				stdout=subprocess.PIPE,
				stderr=subprocess.PIPE)
			
				print flag

				p2 = subprocess.Popen(['samtools', 'bam2fq', '-'],
				stdin=p1.stdout,
				stdout=openfq,
#				print p2
				stderr=subprocess.PIPE)
	#			view_fq1.stdout.close()

				(stdoutdata, stderrdata) = p2.communicate()
				if stdoutdata:
					print 'view_fqrc ', p1.returncode
					print 'bamtofqrc ', p2.returncode
					print "SUCCESS ",stdoutdata
				if stderrdata:
					print 'view_fqrc ', p1.returncode
                        	        print 'bamtofqrc ', p2.returncode
                              		print "ERROR ",stderrdata.strip()

			except OSError as e:
    				print "OSError > ",e.errno
    				print "OSError > ",e.strerror
    				print "OSError > ",e.filename
			except:
    				print "Error > ",sys.exc_info()[0]	
	
	# Count the number of reads in the human filtered fastqs

	with open(filt_fq1) as f1:
               	for i, l in enumerate(f1):
                       	pass

	filt_fw_reads = (i + 1)/4

	print 'Filtered forward reads', filt_fw_reads

        with open(filt_fq2) as f2:
               	for i, l in enumerate(f2):
                       	pass

	filt_rev_reads = (i + 1)/4

	print 'Filtered reverse reads', filt_rev_reads

	result = ET.Element('results')
	items = ET.SubElement(result,'trimming')
	item1 = ET.SubElement(items, 'fw_reads')
	item2 = ET.SubElement(items, 'rev_reads')
	item3 = ET.SubElement(items, 'filt_fw_reads')
	item4 = ET.SubElement(items, 'filt_rev_reads')
	item1.set('result','fwd trimmed reads')
	item2.set('result','rev trimmed reads')
	item3.set('result','fwd filtered reads')
	item4.set('result','rev filtered reads')
	item1.text = str(fw_reads)
	item2.text = str(rev_reads)
	item3.text = str(filt_fw_reads)
	item4.text = str(filt_rev_reads)
	
	xmldata = ET.tostring(result)
	samplexml = '%s_trim_filter.xml' % sample
	xmlfile = open(samplexml, 'w')
	xmlfile.write(xmldata)
	
def split_pops(sample, filt_fq1, filt_fq2):
#	###SPLIT POPULATIONS###
	
	#Essential that this directory is completely separate as splitpops cleans up the directory and deletes all input files
	split_dir = 'splitpops/'
       	os.makedirs(split_dir)

	splitpops = '/home/kieren/Snork6/Snork-PHE/0.6/src/snork.py'
	target_ref = '/phengs/hpc_software/ucl_assembly/new_hcvrefset'
	filt_bam = '%s_filtered.bam' % (sample)

#	#Convert filtered fastqs into bam file
#	filt_bam = '%s_filtered.bam' % (sample)
	try:
		p1 = subprocess.Popen(['java', '-jar',
			       '/phengs/hpc_software/picard-tools/1.111/FastqToSam.jar',
			       'F1='+ filt_fq1,
			       'F2='+ filt_fq2,
		       	       'O='+ filt_bam,
		               'SM='+ sample],	
		stdout=subprocess.PIPE, 
		stderr=subprocess.PIPE)

#	print p1

        	(stdoutdata, stderrdata) = p1.communicate()
		if stdoutdata:
			print "Fastq2bamrc ",p1.returncode
			print "SUCCESS ", stdoutdata
		
		if stderrdata:
			print "Fastq2bamrc ",p1.returncode
        		print "ERROR ",stderrdata.strip()

	except OSError as e:
    		print "OSError > ",e.errno
    		print "OSError > ",e.strerror
    		print "OSError > ",e.filename
	except:
    		print "Error > ",sys.exc_info()[0]

	#Run splitpops.py on bam file
#	/home/kieren/Snork6/Snork-PHE/0.6/src/snork.py splitpops -bin snorktest/src6 -profile None -config /phengs/hpc_software/ucl_assembly/snork.config.wtchg -orgid Hepc -dataid 170926_18 -samplename 170926_18 -bampath temp/170926_18_filtered.bam -targetrefid wtchgR00000071 -targetrefpath /phengs/hpc_software/ucl_assembly/wtchgR00000071.fa -outdir temp/ -logdir temp/ -overwrite False -deleteints True -verbosity DEBUG

#	p2 = subprocess.list2cmdline([splitpops, 'splitpops', '-bin snorktest/src6',
	try:
		p2 = subprocess.Popen([splitpops, 'splitpops',
				     '-bin', 'snorktest/src6',
				     '-profile', 'None',
			             '-config', '/phengs/hpc_software/ucl_assembly/snork.config.wtchg',
			             '-orgid', 'Hepc',
			             '-dataid', sample,
			             '-samplename', sample,
			             '-bampath', filt_bam,
			             '-targetrefid', 'wtchgR00000071',
			             '-targetrefpath', target_ref,
			             '-outdir', split_dir,
			             '-logdir', split_dir,
			             '-overwrite', 'False',
			             '-deleteints', 'True',
			             '-verbosity', 'DEBUG'],

		stdout=subprocess.PIPE,
        	stderr=subprocess.PIPE)
#	print p2
        	(stdoutdata, stderrdata) = p2.communicate()
        	if stdoutdata:
			print "splitpopsrc ",p2.returncode
			print "SUCCESS ",stdoutdata
		if stderrdata:
			print "splitpopsrc ",p2.returncode
			print "ERROR ",stderrdata.strip()
	
	except OSError as e:
    		print "OSError > ",e.errno
    		print "OSError > ",e.strerror
    		print "OSError > ",e.filename
	except:
    		print "Error > ",sys.exc_info()[0]

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

#        	for row in csv.reader(lines, delimiter='\t'):
#                	row[2] = float(row[2])
#                	d[row[0]] = row[2]
#        	for k,v in d.items():
#			print k, '=', v

	result = ET.Element('results')
        items = ET.SubElement(result,'genotypes')
        item1 = ET.SubElement(items, 'genotype1')
        item2 = ET.SubElement(items, 'genotype2')
        item1.set('result',main_geno_name)
        item2.set('result',second_geno_name)
        item1.text = str(main_geno_perc)
        item2.text = str(second_geno_perc)

        xmldata = ET.tostring(result)
        samplexml = '%s_splitpops.xml' % sample
        xmlfile = open(samplexml, 'w')
        xmlfile.write(xmldata)

#	doc_xml = minidom.Document()
        # root of documents for each pipeline module is "result"
#        result = doc_xml.createElement('result')
#        result.setAttribute('sample_id', sample_id)
#        result.setAttribute('type', 'trimmed forward reads')
#        result.setAttribute('value'), fw_reads)



#	#Run bamtools split on splitpops bam
#	bamtools = 'bamtools split -in "temp/"$sample"_filtered.targetpop.readclass.bam" -tag RG:Z ;'
#	os.system(bamtools)
#	#Select largest file from using *RG* tag
#	largest_bam = '`find temp/*RG* -type f | xargs ls -1S | head -n 1`' #Need to check
#	os.system(largest_bam)
#	#Read file for largest bam file filename
#	#largest_bam= cat "temp/"$sample"_splitpop_largest.txt" ;
#
#	#Sort largest bam by readname and then convert selected largest splitpop bam to fastqs
#	sort_largest_bam = 'samtools sort -n $largest_bam_path $largest_bam_path"_sorted"'
#	os.system(sort_largest_bam)
#	bam_to_fastq = 'bamToFastq -i $largest_bam_path"_sorted.bam" -fq "temp/"$sample"_splitpop_1.fastq" -fq2 "temp/"$sample"_splitpop_2.fastq"'
#

def denovo_assembly(sample):
	"""De novo assembly using VICUNA requires the filtered fastq files following human and non-HCV sequence removal"""
	
	###Starting VICUNA###
	#Prepare vicuna config file for specific sample
	
	#Run VICUNA on sample specific config file	
	vicuna_config = '/phengs/hpc_software/ucl_assembly/vicuna_config.txt'
	with open(vicuna_config) as vc:
		cwd = os.getcwd()
		s = vc.read().replace('*DIR*', cwd).replace('*SAMPLE*', sample)
	
	vicuna_config_sample = 'vicuna_config_%s.txt' % (sample)
	with open(vicuna_config_sample, 'w') as vc:
		vc.write(s)

	#Create directory for VICUNA assembly and check to see if it already exists
#	vicuna_dir_path = '%s/assembly/' % (sample)
#	vicuna_dir = os.path.dirname(vicuna_dir_path)
#	if not os.path.exists(vicuna_dir):
#	       	os.makedirs(vicuna_dir)
	
	#Run VICUNA with the new sample config file
	p1 = subprocess.Popen(['vicuna', vicuna_config_sample],
                 stdout=subprocess.PIPE,
                 stderr=subprocess.PIPE)

	(stdoutdata, stderrdata) = p1.communicate()
	print stderrdata
	print stdoutdata
        print 'vicunarc ', p1.returncode

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
	
#	doc_xml = minidom.Document()
        # root of documents for each pipeline module is "result"
#        result = doc_xml.createElement('result')
#        result.setAttribute('sample', sample)

	result = ET.Element('results')

	if int(count) == 0:
		print 'No contig'
		items = ET.SubElement(result, 'Contigs')
		item1 = ET.SubElement(items, 'number of contigs')
		item1.set('Contig value','number')
		item1.text = 'No contigs'    
	

#		result.setAttribute('contigs', 'contigs')
 #       	result.setAttribute('Number', 'No contig')
#		txt = doc_xml.createTextNode('No contig')
#	        result.appendChild(txt)
 #       	doc_xml.appendChild(result)
#		print doc_xml.toprettyxml()

	elif int(count) == 1:
		print 'Single contig'
		items = ET.SubElement(result, 'Contigs')
                item2 = ET.SubElement(items, 'number of contigs')
                item2.set('Contig value','number')
                item2.text = 'Single contig'



#		result.setAttribute('contigs', 'contigs')
#	        result.setAttribute('Number', 'Single contig')
#		txt = doc_xml.createTextNode('Single contig')
 #               result.appendChild(txt)
  #              doc_xml.appendChild(result)
#		print doc_xml.toprettyxml()

	else:
		print 'Number of contigs =', count
		items = ET.SubElement(result, 'Contigs')
                item3 = ET.SubElement(items, 'number of contigs')
                item3.set('Contig value','number')
                item3.text = str(count)




#		result.setAttribute('contigs', 'contigs')
#	        result.setAttribute('Number', str(count))
#		txt = doc_xml.createTextNode(str(count))
 #               result.appendChild(txt)
  #              doc_xml.appendChild(result)
#		print doc_xml.toprettyxml()
		

	print '\n'.join(newdata)
	
	xmldata = ET.tostring(result)
        samplexml = '%s_denovo.xml' % sample
        xmlfile = open(samplexml, 'w')
        xmlfile.write(xmldata) 	
	
	#Very important that this file is named correctly as lastz_analyser.WITH_REVCOMP.pl script will fail later in the pipeline.
	# It must be named samplename.contigs.fasta
 	contigs_increment = 'assembly/%s.contigs.fasta' % (sample)
	
	with open(contigs_increment, 'w') as ci:
		ci.write('\n'.join(newdata))
			
def find_best_ref(sample, hcvfasta):
	"""Following completion of de novo assembly, the best reference sequence from the HCV database needs to
	be selected for mapping to create a draft assembly"""
	
	#Path to current available lastz executable
	lastz_path = '/home/kieren/lastz-distrib/bin/lastz'
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
#		p1.stdout.close()

	(stdoutdata, stderrdata) = p1.communicate()
	print stderrdata
	print stdoutdata
        print 'lastzrc ', p1.returncode

	lastz_log = 'assembly/lastz_besthit.log'
	lastz_bestref = '/phengs/hpc_software/ucl_assembly/lastz_bestref.pl'
	best_ref_fasta = '%s-ref.fasta' % (sample)

	#Run lastz_bestref.pl to find the best matching reference sequence. The perl script requires string concatenation between the flags
	# and arguments (+) for these to run
	p1 = subprocess.Popen(['perl', '-s', lastz_bestref, 
				'-contig_lastz='+ contigs_lastz,
				'-blastdb='+ hcvfasta,
				'-best_ref_fasta='+ best_ref_fasta,
				'-lastz_best_hit_log='+ lastz_log],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE)
        
	(stdoutdata, stderrdata) = p1.communicate()
        print 'BEST REF STDERR', stderrdata
        print 'BEST REF STDOUT', stdoutdata
        print 'lastz_best_refrc ', p1.returncode


def assemble_draft(sample, lastz_path, best_ref_fasta):
	"""The best reference has been defined, now assemble the draft genome"""

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
        print 'Contigs vs bestref STDERR', stderrdata
        print 'Contigs vs bestref STDOUT', stdoutdata
        print 'Contigs_vs_best_refrc ', p1.returncode


	lastz_analyser_rev = '/phengs/hpc_software/ucl_assembly/lastz_analyser.WITH_REVCOMP.pl'
	lastz_analysed_file = 'assembly/lastz_analysed_file'
	lastz_analysed_log = 'assembly/lastz_analyser.log'

	p2 = subprocess.Popen(['perl', '-w', '-s', lastz_analyser_rev,
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
        print 'lastz analyser STDERR', stderrdata
        print 'lastz analyser STDOUT', stdoutdata
        print 'lastz_analyserrc ', p2.returncode

def contig_map(sample, contigs, filt_fq1, filt_fq2, best_ref_fasta, lastz_analysed_file):
	## map reads using BWA MEM to contigs to work out which sequences to choose when contigs overlap

	contigs_sam = 'assembly/%s.contigs.sam' % (sample)
	contigs_bam = 'assembly/%s.contigs.bam' % (sample)
	contigs_sorted = 'assembly/%s.contigs.sorted' % (sample)
	contigs_sorted_final = 'assembly/%s.contigs.sorted.bam' % (sample)
	contigs_mpileup = 'assembly/%s.contigs.mpileup' % (sample)
	genome_maker = '/phengs/hpc_software/ucl_assembly/genome_maker2b.pl'
	genome_fasta = 'assembly/%s-genome.fasta' % (sample)
	genome_maker_log = 'assembly/%s_genome_maker.log' % (sample)
	genome_sam = 'assembly/%s-genome.sam' % (sample)
	genome_bam = 'assembly/%s-genome.bam' % (sample)
	genome_sorted = 'assembly/%s-genome.sorted' % (sample)
	genome_sorted_final = 'assembly/%s-genome.sorted.bam' % (sample)
	genome_mpileup = 'assembly/%s-genome.mpileup' % (sample)
	cons_mv = '/phengs/hpc_software/ucl_assembly/cons_mv.pl'
	cons_pren = 'assembly/%s.consensus1.preNcut.fasta' % (sample)
	genome_mv = 'assembly/%s-genome.fasta.mv' % (sample)
	genome_basefreq = 'assembly/%s-genome.fasta.basefreqs.tsv' % (sample)
	n_remover = '/phengs/hpc_software/ucl_assembly/N_remover_from_consensus.pl'
	cons1_fasta =	'assembly/%s.consensus1.fasta' % (sample)
	cons1_sam = 'assembly/%s.consensus1.sam' % (sample)
	cons1_bam = 'assembly/%s.consensus1.bam' % (sample)
	cons1_sorted = 'assembly/%s.consensus1.sorted' % (sample)
	cons1_sorted_final = 'assembly/%s.consensus1.sorted.bam' % (sample)
	cons1_mpileup = 'assembly/%s.consensus1.mpileup' % (sample)
	cons2_pren = 'assembly/%s.consensus2.preNcut.fasta' % (sample)
	cons1_mv = 'assembly/%s.consensus1.fasta.mv' % (sample)
	cons1_basefreq = 'assembly/%s.consensus1.fasta.basefreqs.tsv' % (sample)
	cons2_fasta = 'assembly/%s.consensus2.fasta' % (sample)
	majvar = '/phengs/hpc_software/ucl_assembly/majvarcheck2_bwa.pl'

	cmd_gm = ['perl', '-w', '-s', genome_maker,
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
 	        print 'bwa index STDERR', stderrdata
        	print 'bwa index STDOUT', stdoutdata
        	print 'bwarc ', p1.returncode

		with open(sam, 'w') as s:
			p1 = subprocess.Popen(['bwa', 'mem',
						'-t', '8', 
						fasta, 
						filt_fq1, 
						filt_fq2],
			stdout=s,
       			stderr=subprocess.PIPE)
	
			(stdoutdata, stderrdata) = p1.communicate()
			print 'bwa mem STDERR', stderrdata
	     	      # print 'bwa mem STDOUT', stdoutdata
        		print 'bwamemrc ', p1.returncode

			p2 = subprocess.Popen(['samtools', 'view', '-bhS',
					'-o', bam, sam],
			stdin=p1.stdout,	
			stdout=subprocess.PIPE,
			stderr=subprocess.PIPE)
	
			(stdoutdata, stderrdata) = p2.communicate()
			print 'samtools view STDERR', stderrdata
        		print 'samtools view STDOUT', stdoutdata
			print 'samtoolsviewrc ', p2.returncode	

	
			p3 = subprocess.Popen(['samtools', 'sort',
        	                                '-@', '8', bam, sorted],
			stdout=subprocess.PIPE,
        		stderr=subprocess.PIPE)

			(stdoutdata, stderrdata) = p3.communicate()
        		print 'samtools sort STDERR', stderrdata
			print 'samtools sort STDOUT', stdoutdata
        		print 'samtoolssortrc ', p3.returncode

		
			p4 = subprocess.Popen(['samtools', 'index', final],
        	        stdout=subprocess.PIPE,
                	stderr=subprocess.PIPE)

  	                (stdoutdata, stderrdata) = p4.communicate()
        	        print 'samtools index STDERR', stderrdata
                	print 'samtools index STDOUT', stdoutdata
               	        print 'samtoolsindexrc ', p4.returncode
		
			p5 = subprocess.Popen(['samtools', 'mpileup', '-f',
						fasta,
						'-d', '1000000',
						'-o', mpileup, final],
	                stdout=subprocess.PIPE,
        	        stderr=subprocess.PIPE)

                	(stdoutdata, stderrdata) = p5.communicate()
	                print 'samtools mpileup STDERR', stderrdata
#       	        print 'samtools mpileup STDOUT', stdoutdata
	                print 'samtoolsmpileuprc ', p5.returncode
	
		#Helpful method for checking exactly what is sent to the command line (subprocess.list2cmdline)
			if fasta==contigs:	
			#	p6 = subprocess.list2cmdline(cmd_gm)
				p6 = subprocess.Popen(cmd_gm, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	
		#Print the full command line argument
				print 'fasta is contigs'

				(stdoutdata, stderrdata) = p6.communicate()
	        		print 'gm or cons STDERR', stderrdata
				print 'gm or cons STDOUT', stdoutdata
	       			print 'gm or cons rc ', p6.returncode

			elif fasta==genome_fasta:

				p7 = subprocess.Popen( ['perl', '-w', '-s', cons_mv,
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
        	                print 'Consmv STDERR', stderrdata
                	        print 'Consmv STDOUT', stdoutdata
                        	print 'Consmv rc ', p7.returncode
				
				p8 = subprocess.Popen(['perl', '-w', '-s', n_remover,
                                     			'-cutoff=46',
                                     			cons_pren],

				stdout=open(cons1_fasta, 'w'), 
				stderr=subprocess.PIPE)				
				
                                (stdoutdata, stderrdata) = p8.communicate()
                                print 'N_remover STDERR', stderrdata
                                print 'N_remover STDOUT', stdoutdata
                                print 'N_remover rc ', p8.returncode

			else:

				p9 = subprocess.Popen( ['perl', '-w', '-s', cons_mv,
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
                                print 'Consmv STDERR', stderrdata
                                print 'Consmv STDOUT', stdoutdata
                                print 'Consmv rc ', p9.returncode

                                p10 = subprocess.Popen(['perl', '-w', '-s', n_remover,
                                                        '-cutoff=46',
                                                        cons2_pren],

                                stdout=open(cons2_fasta, 'w'),
                                stderr=subprocess.PIPE)

                                (stdoutdata, stderrdata) = p10.communicate()
                                print 'N_remover STDERR', stderrdata
                                print 'N_remover STDOUT', stdoutdata
                                print 'N_remover rc ', p10.returncode

				p11 = subprocess.Popen(['perl', '-w', '-s', majvar,
                                       		       '-mvpath='+ cons1_mv,
                                        	       '-basefreq='+ cons1_basefreq,
                                        	       '-fwdreads='+ filt_fq1,
		                                       '-revreads='+ filt_fq2], 
				stdout=subprocess.PIPE, 
				stderr=subprocess.PIPE)
				
				(stdoutdata, stderrdata) = p11.communicate()
 	                        print 'Majvar STDERR', stderrdata
  	                        print 'Majvar STDOUT', stdoutdata
   	                        print 'Majvar rc ', p11.returncode
				
				
				p12 = subprocess.Popen(['samtools', 'view',
						       '-S',
                                                       '-F0x4',
                                                       '-c', contigs_sam],
				stdout=subprocess.PIPE,
                                stderr=subprocess.PIPE)

                                (stdoutdata, stderrdata) = p12.communicate()
                                print 'samtools view mapcount STDERR', stderrdata
                                print 'samtools view mapcount STDOUT', stdoutdata
                                print 'samtools view mapcount rc ', p12.returncode

				result = ET.Element('results')
        			items = ET.SubElement(result,'mapping')
        			item1 = ET.SubElement(items, 'mapping')
        			item1.set('result','reads mapped to contigs')
	        		item1.text = str(stdoutdata)

        			xmldata = ET.tostring(result)
        			samplexml = '%s_mapping.xml' % sample
        			xmlfile = open(samplexml, 'w')
        			xmlfile.write(xmldata)

def quasibam(sample, cons1_sorted_final):
#	"""Run quasibam on consensus bam file and consensus fasta sequence"""

	consensus = 'assembly/%s.consensus2.fasta' % (sample)	
#	quasifas = '%s/assembly/%s_vicuna_bwa_quasibam.fas' % (sample, sample)
#	quasitxt = '%s/assembly/%s_vicuna_bwa_quasibam.txt' % (sample, sample)
#	quasierr = '%s/assembly/%s_vicuna_bwa_quasibam.err' % (sample, sample)

#	with open(quasifas, 'w') as qf:

#	p1 = subprocess.list2cmdline(['quasi_bam',
#                                       cons1_sorted_final,
#                                       consensus])
#	print p1

	p1 = subprocess.Popen(['quasi_bam',
				cons1_sorted_final, 
				consensus,
				'-c', '15'],
#				'-o1', qf],
#				'-o2', quasitxt,
#				'-o3', quasierr],
                
	stdout=subprocess.PIPE,
	stderr=subprocess.PIPE)

	(stdoutdata, stderrdata) = p1.communicate()
      	print stderrdata
	print stdoutdata
      	print 'quasibamrc ', p1.returncode


if __name__ == '__main__':

        main()

