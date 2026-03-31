#!/usr/bin/python

import sys
import os
import pandas
import re
from Bio import SeqIO
import gzip

# Arguments
print('Number of arguments: ', len(sys.argv), 'arguments.')
print('Input read directory: ', sys.argv[1])
print('Statistics file: ', sys.argv[2])
print('Sample: ', sys.argv[3])
print('\n---\n')

read_dir = sys.argv[1]
stats_file =  sys.argv[2]
sampleID = sys.argv[3]

forward_fastq = os.path.join(read_dir, sampleID + '_R1.clean.EC.fastq.gz')
reverse_fastq = os.path.join(read_dir, sampleID + '_R2.clean.EC.fastq.gz')

F_count = 0
F_length_sum = 0
F_GC_sum = 0
F_GC_perc_sum = 0
R_count = 0
R_length_sum = 0
R_GC_sum = 0
R_GC_perc_sum = 0

with gzip.open(forward_fastq, 'rt') as input_fastq:
	read_seq = SeqIO.parse(input_fastq, 'fastq')
	for record in read_seq:
		F_count = F_count + 1
		seq_length = len(str(record.seq))
		F_length_sum = F_length_sum + seq_length
		GC_count = str(record.seq).count('G') + str(record.seq).count('C') + str(record.seq).count('g') + str(record.seq).count('c')
		F_GC_sum = F_GC_sum + GC_count
		GC_perc = GC_count/seq_length*100
		F_GC_perc_sum = F_GC_perc_sum + GC_perc
F_GC_perc = F_GC_perc_sum/F_count
print(str(F_count) + ' forward reads at ' + str(F_GC_perc) + '% average GC')

with gzip.open(reverse_fastq, 'rt') as input_fastq:
       	read_seq = SeqIO.parse(input_fastq, 'fastq')
       	for record in read_seq:
                R_count	= R_count + 1
               	seq_length = len(str(record.seq))
               	R_length_sum = R_length_sum + seq_length
               	GC_count = str(record.seq).count('G') + str(record.seq).count('C') + str(record.seq).count('g') + str(record.seq).count('c')
                R_GC_sum = R_GC_sum + GC_count
                GC_perc	= GC_count/seq_length*100
                R_GC_perc_sum =	R_GC_perc_sum +	GC_perc
R_GC_perc = R_GC_perc_sum/R_count
print(str(R_count) + ' reverse reads at ' + str(R_GC_perc) + '% average GC')

with open(stats_file, 'a') as outfile:
	outfile.write('\t'.join([sampleID, 
                                 str(F_count), str(F_length_sum), str(F_GC_sum), str(F_GC_perc),
                                 str(R_count), str(R_length_sum), str(R_GC_sum), str(R_GC_perc)]) + '\n')
