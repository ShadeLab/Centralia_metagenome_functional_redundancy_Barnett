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
contig_cov = pd.DataFrame()
sub_annotation = pd.read_csv(os.path.join(cov_dir, sample_list[0] + '.cov.txt'), sep='\t')
sub_annotation['MAG_num'] = sub_annotation['#ID'].str.split('_').str[0]
contig_cov['MAG_num'] = list(set(sub_annotation['MAG_num']))

for sampleID in sample_list:
	sub_annotation = pd.read_csv(os.path.join(cov_dir, sampleID + '.cov.txt'), sep='\t')
	sub_annotation['MAG_num'] = sub_annotation['#ID'].str.split('_').str[0]
	sub_annotation[sampleID] = sub_annotation['Plus_reads'] + sub_annotation['Minus_reads']
	sub_annotation = sub_annotation.groupby(['MAG_num'])[sampleID].sum()
	contig_cov = contig_cov.merge(sub_annotation, how='left', on=['MAG_num'])
	sub_annotation = None


############################################################
###################### Writing output ######################
############################################################

contig_cov.to_csv(output_file, header=True, index=False, sep='\t')
