#!/usr/bin/env python3
# module load hmmer/latest

import os
import argparse
import numpy as np
import pandas as pd
from subprocess import Popen
from Bio.SeqIO.FastaIO import SimpleFastaParser as sfp
from itertools import islice

hmm = '/phengs/hpc_storage/home/dbibby/hmmer/HCV/HCV_domains'
DEVNULL = open(os.devnull, 'wb')

def main(tabular, fas, threshold, domains, lengths, percentiles):
	"""
	Uses HMMER to query a Whole Genome Sequence of HCV agains a set of HMM
	profiles to obtain co-ordinates for domains of HCV (defaults to NS3,
	NS5a and NS5b only, although 5'UTR, Core, E1, E2, and NS4b are possible).
	Resulting co-ordinates are used to derive statistics of depth and coverage
	from a quasibam nucleotide frequency table as well as the domain sequences
	themselves.
	For resistance purposes, the domains can be 3'-limited, ie. a terminal
	amino acid to be included can be specified. This is useful in e.g. NS5a,
	where only the first 213 amino acids are analysed.
	
	arguments
	---------
	qb:      file path to the quasibam nucleotide frequency file
	fas:     file path to the reference sequence used to generate <qb>
	domains: tuple containing text labels for the required domains.
	lengths: tuple of same length as domains containing the last amino acid
	         positions to be counted in the stats (eg. NS5a to 213)
	         If None, the entire domains are included.
	
	outputs
	-------
	
	"""
	def stats(hmmer_line, length):
	    #   where domain is missing, abort analysis
		if not hmmer_line: return None
		#	locates the co-ordinates of the domain
		start, end = map(int, hmmer_line[10:12])
		#	defaults to entire domain unless specified
		if not length: length = 1 + end - start
		
		#	Reference nucleotide and depth are conveniently located in
		#	adjacent columns of the QuasiBAM table (i.e. [1:3])
		ref, depth = islice(zip(*quasi_table[start:][:length]), 1, 3)

		#	Turns the list into a sequence string
		ref = ''.join(ref)
		
		#	Creates a numpy array of the depths
		depth = np.fromiter(map(int, depth),
							dtype=np.int64)
		
		#	Calculates the proportion of sites having a depth greater than the
		#	specified threshold (default=100)
		coverage = round(100 * depth[depth >= threshold].size / depth.size, 2)
		
		#	Generates an argsort array of depth to enable easy calculation of
		#	percentiles
		as_d = np.argsort(depth, kind='mergesort')		
		depth_q = list(map(lambda x: depth[as_d[int(depth.size * x / 100)]],
					   percentiles))
		
		return [start, end, ref, coverage] + depth_q
	
	#	Loads the QuasiBAM table into memory
	quasi_table = [i.split() for i in open(tabular).readlines()]
	#   Initialises a pandas dataframe to receive the output data
	df = pd.DataFrame(columns=['Start', 'End', 'Ref', 'Coverage'] + percentiles,
					  index=domains)
					  
	#	Uses nhmmscan in HMMER to align the input fasta to a collection of
	#	HMM profiles stored in 'hmm'. These have been generated from the
	#	LANL curated alignments, delimited by each domain.
	#	stdout is suppressed as the summary table produced by --dfamtblout is
	#	sufficient
	hmmer_out = 'table.out'
	HMMER_call = Popen(['nhmmscan',
						'--dfamtblout', hmmer_out,
						hmm, fas],
		stdout=DEVNULL)
	HMMER_call.wait()
	
	
	with open(hmmer_out) as result:
		#	loads the output lines in to a dictionary to enable the order of
		#	input domains & lengths to be preserved
		h_dict = {}
		for i in map(str.split, result):
			h_dict[i[0]] = h_dict.get(i[0], i[1:])
		
		#	obtains the stats for each domain
		for d, l in zip(domains, lengths):
		    #   writes the output of stats to the domain row of the dataframe
			df.loc[d] = stats(h_dict.get(d, None), l * 3)
	
	#   writes the dataframe as a csv
	df.to_csv(os.path.splitext(fas)[0] + '_domains.csv')

if __name__ == '__main__':
	parser = argparse.ArgumentParser(prog='HCV_domains.py')
	parser.add_argument('-q', metavar='tabular file',
						required=True,
						help='quasibam nucleotide frequency table')
	
	parser.add_argument('-f', metavar='FASTA',
						required=True,
						help='fasta file of reference sequence')
	
	parser.add_argument('-d', metavar='domains',
						required=False,
						default='NS3,NS5a,NS5b',
						help='Comma-separated list of HCV domains to query. ' 
						'Defaults to NS3, NS5a, NS5b')
	
	parser.add_argument('-l', metavar='lengths',
						required=False,
						default='0,0,0',
						help='Comma-separated list of terminal amino acid '
						'positions to consider in the domains specified by -d')
	
	parser.add_argument('-t', metavar='threshold',
						required=False,
						type=int,
						default=100,
						help='Depth threshold for inclusion in coverage stat')
	
	parser.add_argument('-p', metavar='percentiles',
						required=False,
						default='25,50,75',
						help='Comma-separated list of quartile percentages')
	
	args = parser.parse_args()
	q = args.q
	f = args.f
	t = args.t
	d = args.d.split(',')
	l = list(map(int, args.l.split(',')))
	p = list(map(int, args.p.split(',')))
	
	main(tabular=q,
		 fas=f,
		 threshold=t,
		 domains=d,
		 lengths=l,
		 percentiles=p)
		
