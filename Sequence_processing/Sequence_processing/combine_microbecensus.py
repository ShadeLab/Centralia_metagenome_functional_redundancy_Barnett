#!/usr/bin/python

import sys
import os
import re

with open('./sample_list.txt', 'r') as sam_file:
	sam_list = sam_file.readlines()
	sam_list = [line.rstrip() for line in sam_list]

with open('/mnt/research/ShadeLab/Barnett/Centralia_metagenome/annotations/microbecensus/microbecensus_comb_out.txt', 'w') as outfile:
	outfile.write('SequenceID\taverage_genome_size\ttotal_bases\tgenome_equivalents\n')
	for samID in sam_list:
		mc_file = '/mnt/research/ShadeLab/Barnett/Centralia_metagenome/annotations/microbecensus/' + samID + '_census.txt'
		with open(mc_file, 'r') as infile:
			lines = infile.readlines()
			lines = [line.rstrip() for line in lines]
		outfile.write(samID + '\t')
		for line in lines:
			if line.startswith('average_genome_size'):
				outfile.write(line.split(':\t')[1] + '\t')
			elif line.startswith('total_bases'):
				outfile.write(line.split(':\t')[1] + '\t')
			elif line.startswith('genome_equivalents'):
				outfile.write(line.split(':\t')[1])
		outfile.write('\n')
