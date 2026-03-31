#!/usr/bin/python

import sys
import os
import re
from Bio import SeqIO
import pandas as pd

# Arguments
print('Number of arguments: ', len(sys.argv), 'arguments.')
print('MAG list file: ', sys.argv[1])
print('MAG fasta directory: ', sys.argv[2])
print('Combined fasta file: ', sys.argv[3])
print('MAG mapping file: ', sys.argv[4])
print('\n---\n')

MAG_list_file = sys.argv[1]
MAG_fasta_dir =  sys.argv[2]
comb_fasta_file = sys.argv[3]
mapping_file = sys.argv[4]

MAG_df = pd.read_csv(MAG_list_file, sep=',')
MAG_list = list(MAG_df['genome'])

MAG_num = 0
with open(comb_fasta_file, 'w') as out_fasta:
	with open(mapping_file, 'w') as out_map:
		out_map.write('MAG_ID\tContigID\tMapped_ID\n')
		for MAG in MAG_list:
			MAG_num += 1
			contig_num = 0
			MAG_fasta = SeqIO.parse(os.path.join(MAG_fasta_dir, MAG), 'fasta')
			for record in MAG_fasta:
				contig_num += 1
				out_fasta.write('>MAG' + str(MAG_num) + '_contig' + str(contig_num) + '\n') 
				out_fasta.write(str(record.seq) + '\n')
				out_map.write(MAG + '\t' + record.id + '\t' + 'MAG' + str(MAG_num) + '_contig' + str(contig_num) + '\n')
