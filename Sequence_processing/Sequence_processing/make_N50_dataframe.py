#!/usr/bin/python

import sys
import os
import pandas
import re

with open('./S6_calculate_N50.out', 'r') as infile:
	lines = infile.readlines()
	lines = [line.rstrip() for line in lines]

with open('./N50_data.txt', 'w') as outfile:
	outfile.write('SequenceID\tN50\tContig_count\tTotal_length\tLongest_contig\n')
	for line in lines:
		if line.startswith('Cen'):
			outfile.write(line + '\t')
		elif line.startswith('N50'):
			outfile.write(line.split(': ')[1] + '\t')
		elif line.startswith('Sequences'):
                        outfile.write(line.split(': ')[1] + '\t')
		elif line.startswith('Total length'):
                        outfile.write(line.split(': ')[1] + '\t')
		elif line.startswith('Longest sequence'):
                        outfile.write(line.split(': ')[1] + '\t')
		elif line.startswith('Done!'):
			outfile.write('\n')
