"""
File:		find_seq_in_contigs.py
Aythor:		Somayeh Kafaie
Date:		September 2020
Purpose:	To find the corresponding contigs for sequences

To run:
	python find_seq_in_contigs.py --contig/-C <contig file path>
    --seq/-s <the sequence file>

"""

################################################################################

import sys
import os
import gfapy
import re
import argparse
import difflib
import datetime
import csv
import collections
from Bio import SeqIO

OUTPUT_DIR = 'temp'
result_file = OUTPUT_DIR+'/sequence_match_list'+datetime.datetime.now().strftime('%Y-%m-%d_%H-%M')+'.csv'

def write_results(result_file, query, contig_list):
	"""
	To write the list of contigs found in the assembly graph (contig_list) for
	a given sequence (query) in a file (result_file)
	Parameters:
	 	result_file: The file to be written in
		query: the sequence
		contig_list: the list of contigs
	"""
	with open(result_file, mode='a+') as file:
		common_writer = csv.writer(file, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
		if (contig_list):
			for contig in contig_list:
				common_writer.writerow([query, len(query)-1, contig.name, contig.identity, contig.matched_length,\
					contig.q_start, contig.q_end, contig.c_start, contig.c_end])
		else:
			common_writer.writerow([query, len(query)-1, '', '', '', '', '', '', ''])

def find_sequence_match(query, contig_file, out_dir = OUTPUT_DIR, blast_ext = ''):
	"""
	To check if the sequence (query) can be found in the list of contigs
	# contig_file is considered our database and we run blast over it
	Parameters:
		query: the sequence we look for in the contigs
		contig_file: the file containing contigs
		out_dir: the output directory to store all created files there
		blast_ext: an extension used for naming blast file
	Return:
		the list of matched contigs and some further information including identity,
		matched_length, the start and end of matched part of the query, and the start
		and end of matched part of the contig
	"""
	#Create DB from contig file
	command = 'makeblastdb -in '+contig_file +' -parse_seqids -dbtype nucl'
	os.system(command)
	#write the sequence into a fasta file
	query_file = out_dir+'/query.fasta'
	file = open(query_file, 'w')
	file.write(query)
	file.close()
	#run blast query for alignement
	blast_file_name = out_dir+'/blast'+blast_ext+'.csv'
	command = 'blastn -query '+query_file+' -db '+contig_file+\
		' -task blastn -outfmt 10 -max_target_seqs 5 -evalue 0.5 -perc_identity 95 > '+ blast_file_name

	os.system(command)
	contig_list = []
	contig = collections.namedtuple('contig', 'name identity matched_length q_start q_end c_start c_end')
	with open(blast_file_name, 'r') as file1:
		myfile = csv.reader(file1)
		for row in myfile:
			mycontig = contig(name=row[1], identity=row[2], matched_length=row[3], q_start=row[6],
		 					q_end=row[7], c_start=row[8], c_end=row[9])
			contig_list.append(mycontig)

	return contig_list



def find_sequences(seq_file, contig_file):
	"""
	To extract the sequences, find the list of contigs containing each sequence
	and write the results in a file
	Parameters:
		seq_file:	the file containing the sequences we look for in contigs
		contig_file:the file containing contigs
	"""
	#create the output directory
	if not os.path.exists(OUTPUT_DIR):
		os.makedirs(OUTPUT_DIR)
	#write headers into the result file
	with open(result_file, mode='w') as file:
		result_writer = csv.writer(file, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
		result_writer.writerow(['sequence','sequence_length','contig','identity',\
		 			'matched_length', 'seq_start', 'seq_end', 'contig_start', 'contig_end'])
	#extract sequences and align them
	with open(seq_file, 'r') as myfile:
		for line in myfile:
			if line.startswith('>') or line.startswith('Path') or line.startswith('The'):
				continue
			else:
				contig_list = find_sequence_match(line, contig_file)
				write_results(result_file, line, contig_list)


def main(args):
	"""
	"""
	if args.seq and args.contig:
		find_sequences(args.seq, args.contig)
	else:
		print('please enter the path for the sequence file and contig file')
		sys.exit()


if __name__=="__main__":

	text = 'This code is used to find the sequences in the contig file'
	parser = argparse.ArgumentParser(description=text)
	parser.add_argument('--contig', '-C', type=str, default='',
		help='the path of the file containing the sequence of all contigs')
	parser.add_argument('--seq','-S', type=str, default = '',
		help = 'the path of the file containing the sequences')
	args = parser.parse_args()
	main(args)
