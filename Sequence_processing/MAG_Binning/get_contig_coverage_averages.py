#!/usr/bin/python

import sys
import os
import pandas as pd
import re
from Bio import SeqIO

# Arguments
print('Number of arguments: ', len(sys.argv), 'arguments.')
print('Coverage directory: ', sys.argv[1])
print('Sample list: ', sys.argv[2])
print('Output file: ', sys.argv[3])
print('\n---\n')

cov_dir = sys.argv[1]
sample_list_file =  sys.argv[2]
output_file = sys.argv[3]

with open(sample_list_file, 'r') as f:
	sample_list = [line.rstrip() for line in f]


############################################################
#################### Coverage annotations ##################
############################################################

print('Getting coverage outputs')

with open(output_file, 'w') as outfile:
	outfile.write('SequenceID\tAvg_coverage\n')
	for sampleID in sample_list:
		sub_annotation = pd.read_csv(os.path.join(cov_dir, sampleID, 'work_files/metabat_depth.txt'), sep='\t')
		sub_avg_coverage = sum(sub_annotation['totalAvgDepth'])/len(sub_annotation['totalAvgDepth'])
		sub_annotation = None
		outfile.write(sampleID + '\t' + str(sub_avg_coverage) + '\n')
