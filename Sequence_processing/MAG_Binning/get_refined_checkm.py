#!/usr/bin/python

import sys
import os
import pandas as pd
import re
from Bio import SeqIO

# Arguments
print('Number of arguments: ', len(sys.argv), 'arguments.')
print('Sample list: ', sys.argv[1])
print('Refined bin directory: ', sys.argv[2])
print('Combined directory: ', sys.argv[3])
print('\n---\n')

sample_list_file = sys.argv[1]
bin_dir =  sys.argv[2]
comb_dir = sys.argv[3]

with open(sample_list_file, 'r') as f:
	sample_list = [line.rstrip() for line in f]


############################################################
#################### Combine checkM data ###################
############################################################

print('Getting checkM outputs from refined bins')
checkM_df = pd.DataFrame()
for sampleID in sample_list:
	stats_file = os.path.join(bin_dir, sampleID, 'metawrap_50_10_bins.stats')
	if os.path.exists(stats_file):
		sub_annotation = pd.read_csv(stats_file, sep='\t')
		sub_annotation['genome'] = sampleID + '_' + sub_annotation['bin'] + '.fa'
		sub_annotation = sub_annotation[['genome', 'completeness', 'contamination', 'GC', 'lineage', 'N50', 'size', 'binner']]
		checkM_df = checkM_df.append(sub_annotation, ignore_index=True)
		sub_annotation = None
	else:
		print(sampleID + ' is missing its CheckM stats file')

checkM_df.to_csv(os.path.join(comb_dir, 'refined_bin_checkM.txt'), header=True, index=False, sep='\t')
