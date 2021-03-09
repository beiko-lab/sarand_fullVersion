"""
File:		full_pipeline.py
Aythor:		Somayeh Kafaie
Date:		September 2020


Purpose:	To extract the neighborhood of an AMR gene from an assembly graph,
			annotate it and compare with the reference genome(s)

To run:
	python full_pipeline.py

NOTE: all parameters can be set in params.py
NOTE: if use_RGI = TRUE, make sure either RGI has been installed system-wide or
	you already are in the environment RGI installed in!
"""

import sys
import os
import copy
import gfapy
import re
import argparse
import difflib
import datetime
import csv
from csv import DictReader
import collections
from Bio import SeqIO
from gfapy.sequence import rc
import enum
import subprocess
import random
import shutil
from dna_features_viewer import GraphicFeature, GraphicRecord
import matplotlib.pyplot as plt
import numpy as np
from PIL import Image
import json
import logging
from functools import partial
from multiprocessing.pool import Pool
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd

from params import Pipeline_tasks, Insertion_type
from extract_neighborhood import neighborhood_graph_extraction, neighborhood_sequence_extraction,\
							extract_amr_neighborhood_in_ref_genome, extract_nodes_in_path,\
							compare_two_sequences
from find_seq_in_contigs import find_sequence_match
from annotation_visualization import visualize_annotation
from utils import initialize_logger, check_reads, str2bool, print_params, verify_file_existence,\
			retrieve_AMR, extract_files


ASSEMBLY_FILE = 'assembly_graph_with_scaffolds.gfa'
CONTIG_FILE = 'contigs.fasta'
ALL_AMR_SEQUENCES ='nucleotide_fasta_protein_homolog_model_without_efflux.fasta'
AMR_FAMILY_INFO = 'aro_index.tsv'
AMR_DIR_NAME = 'AMR/'
AMR_SEQ_DIR = 'sequences/'
AMR_ALIGN_DIR = 'alignments/'
AMR_OVERLAP_FILE = 'overlaps.txt'
SUBGRAPH_DIR_NAME = 'subgraphs/'
SEQ_DIR_NAME = 'sequences_info'
SEQ_NAME_PREFIX = 'ng_sequences_'
ANNOTATION_DIR = 'annotations'
EVAL_DIR = 'evaluation'

#1 means we annotate the contig file first and treat any gene from the same family of target AMR as AMR itself
#2 means we annotate the contig file first and only look for the exact name of target AMR
#3 means we first align our AMR into contig file and after extracting neighbothood sequences annotate them
CAMI_REF_NEIGHBORHOOD_METHOD = 3

MULTIPLE_SEQ_LENGTH = False
MULTIPLE_COV_THR = True
#METASPADES_USED = True

def validate_task_values(tasks):
	"""
	To check if the task(s) entered by the user are valid
	Parameter:
		tasks: it's either a number representing the task number valid in Pipeline_tasks
			or two numbers denoting the start and end task which are valid values in Pipeline_tasks
	Return:
		the list of tasks to be done
		For example if tasks =[1, 4] return [1, 2, 3, 4]
		special case: tasks = [0] return [1, 2, 3, 4, 5, 6]
	"""
	task_error_message = "For the entire pipeline choose "+str(Pipeline_tasks.all.value)+"; otherwise\
	either provide a number representing one of the following tasks or two numbers\
	to denote the start and end tasks (and of course all tasks in the middle will be run).\n \
	Here is the list:\nmetagenome_creation = "+str(Pipeline_tasks.metagenome_creation.value)+\
	"\nread_simulation = "+str(Pipeline_tasks.read_simulation.value)+\
	"\nassembly = "+str(Pipeline_tasks.assembly.value)+\
	"\ngraph_neighborhood = "+str(Pipeline_tasks.graph_neighborhood.value)+\
	"\nsequence_neighborhood = "+str(Pipeline_tasks.sequence_neighborhood.value)+\
	"\nneighborhood_annotation = "+str(Pipeline_tasks.neighborhood_annotation.value)+\
	"\nneighborhood_evaluation = "+str(Pipeline_tasks.neighborhood_evaluation.value)
	#"\ncontig_matching = "+str(Pipeline_tasks.contig_matching.value)+\
	task_list = []
	if len(tasks) > 2:
		logging.error("ERROR: There are more than two numbers in the task list!\n" + task_error_message)
		#print("ERROR: There are more than two numbers in the task list!\n" + task_error_message)
		sys.exit()

	valid_task_values = [item.value for item in Pipeline_tasks]
	for task in tasks:
		if int(task) not in valid_task_values:
			logging.error("ERROR: invalid task number(s)!\n" + task_error_message)
			#print("ERROR: invalid task number(s)!\n" + task_error_message)
			sys.exit()

	if len(tasks)==2 and int(tasks[0])>int(tasks[1]):
		logging.error("ERROR: The first task number should be smaller than the second task\
		 in the list!\n" + task_error_message)
		#print("ERROR: The first task number should be smaller than the second task\
		 #in the list!\n" + task_error_message)
		sys.exit()

	if len(tasks)==1 and int(tasks[0])==Pipeline_tasks.all.value:
		return valid_task_values
	if len(tasks)==1:
		return [int(tasks[0])]
	for task in list(range(int(tasks[0]), int(tasks[1])+1)):
		task_list.append(task)
	return task_list

def read_amr_alignment(tsv_file, amr_name, threshold =  99):
	paths_info = []
	found = False
	with open(tsv_file) as tsvfile:
		reader = csv.reader(tsvfile, delimiter='\t')
		#skip the header
		next(reader)
		for row in reader:
			coverage = float(re.sub('[%]','',row[3]))
			identity = float(re.sub('[%]','',row[5]))
			if int(coverage) >= threshold and int(identity)>=threshold:
				found = True
				cell_info = row[1].strip()
				nodes, orientation_list, start_pos, end_pos = extract_nodes_in_path(cell_info)
				path_info = {'nodes':nodes, 'orientations':orientation_list,
								'start_pos':start_pos, 'end_pos':end_pos}
				paths_info.append(path_info)
	if not found:
		logging.onfo('ERROR: AMR "'+amr_name+'" has an alignment file but no path-info was found')
	return found, paths_info

def extract_amr_family_info(file_name = AMR_FAMILY_INFO):
	"""
	"""
	family_list = []
	family_info = []
	with open(file_name, 'r') as myfile:
		myreader = DictReader(myfile, delimiter='\t')
		for row in myreader:
			gene = row['Model Name'].strip().replace(' ','_')
			gene_name = ''.join(e for e in gene if e.isalpha() or e.isnumeric() or e=='_' or e=='-')
			if row['AMR Gene Family'] not in family_list:
				family_list.append(row['AMR Gene Family'])
				family_info.append([gene_name.lower()])
			else:
				myindex = family_list.index(row['AMR Gene Family'])
				family_info[myindex].append(gene_name.lower())

	return {'family':family_list, 'gene_list':family_info}

# def write_info_in_annotation_file(annotation_writer, visual_annotation_writer,
# 					seq_description, seq, contig_name, seq_info, coverage, use_RGI, found):
def write_info_in_annotation_file(annotation_writer, visual_annotation_writer,
								gene_info, use_RGI, found, len_seq = None):
	"""
	To write annotation details into files
	Parameters:
		annotation_writer:	annotation file containing all annotations
		visual_annotation_writer: annotation file containing unique annotations
		seq_description: a small description of the sequence used for naming
		seq: the extracted dna sequence that has been annotated
		gene_info: annotation info
		contig_name: the name of contig matching this extracted sequence (if there is any contig)
		use_RGI: True if RGI has been used to annotate AMRs
		found: True if the annotation info has already found in other annotated sequences
	"""
	seq = gene_info['seq_value']
	seq_description = gene_info['seq_name']
	if len_seq is None:
		len_seq = len(seq)
	if use_RGI:
		annotation_writer.writerow([seq_description, seq, len_seq,
							gene_info['gene'], gene_info['prokka_gene_name'],
							gene_info['product'], gene_info['length'],
							gene_info['start_pos'], gene_info['end_pos'],
							gene_info['RGI_prediction_type'], gene_info['coverage'],
							gene_info['family'], gene_info['target_amr']])
		if not found:
			visual_annotation_writer.writerow([seq_description, seq, len_seq,
								gene_info['gene'], gene_info['prokka_gene_name'],
								gene_info['product'], gene_info['length'],
								gene_info['start_pos'], gene_info['end_pos'],
								gene_info['RGI_prediction_type'], gene_info['coverage'],
								gene_info['family'], gene_info['target_amr']])
	else:
		annotation_writer.writerow([seq_description, seq, len_seq,
							gene_info['gene'],
							gene_info['product'], gene_info['length'],
							gene_info['start_pos'], gene_info['end_pos'],
							gene_info['coverage'], gene_info['target_amr']])
		if not found:
			visual_annotation_writer.writerow([seq_description, seq, len_seq,
								gene_info['gene'],
								gene_info['product'], gene_info['length'],
								gene_info['start_pos'], gene_info['end_pos'],
								gene_info['coverage'], gene_info['target_amr']])
	# if use_RGI:
	# 	annotation_writer.writerow([seq_description, seq, len(seq),
	# 						contig_name, gene_info['gene'], gene_info['prokka_gene_name'],
	# 						gene_info['product'], gene_info['length'],
	# 						gene_info['start_pos'], gene_info['end_pos'],
	# 						gene_info['RGI_prediction_type'], coverage])
	# 	if not found:
	# 		visual_annotation_writer.writerow([seq_description, seq, len(seq),
	# 							contig_name, gene_info['gene'], gene_info['prokka_gene_name'],
	# 							gene_info['product'], gene_info['length'],
	# 							gene_info['start_pos'], gene_info['end_pos'],
	# 							gene_info['RGI_prediction_type'], coverage])
	# else:
	# 	annotation_writer.writerow([seq_description, seq, len(seq),
	# 						contig_name, gene_info['gene'],
	# 						gene_info['product'], gene_info['length'],
	# 						gene_info['start_pos'], gene_info['end_pos'], coverage])
	# 	if not found:
	# 		visual_annotation_writer.writerow([seq_description, seq, len(seq),
	# 							contig_name, gene_info['gene'],
	# 							gene_info['product'], gene_info['length'],
	# 							gene_info['start_pos'], gene_info['end_pos'], coverage])


def insert_amr_in_location(genome_file, amr_seq, location, output_dir):
	"""
	To insert The AMR sequence in a given location (attaching to the end of a given line)
	of the genome sequence
	Parameters:
		genome_file: the file in which the AMR file is to be inserted
		amr_seq:	the AMR sequence to be inserted
		location:	the line number to insert the AMR gene after it
		output_dir:	the output directory to store the result
	Return:
		the address of the generated file after the process
	"""
	base=os.path.basename(genome_file)
	new_file_name = output_dir+os.path.splitext(base)[0]+'_'+str(location)+os.path.splitext(base)[1]
	new_file = open(new_file_name, 'w')
	with open(genome_file, 'r') as read_obj:
		for count,line in enumerate(read_obj):
			if count==location:
				line = line[:-1] + amr_seq
			new_file.write(line)
	new_file.close()
	return new_file_name

def create_metagenome_with_amr_insertion(ref_genome_files, number_of_insertions,insertion_type,
						insertion_locations, amr_file, output_dir, metagenome_file = 'metagenome.fasta'):
	"""
	To create a metagenome sample from some genomes after inserting the AMR gene in them
	Parameters:
		ref_genome_files:	the list of reference genomes in which AMR gene is supposed to be inserted
		number_of_insertions:	the number of times the AMR gene is supposed o be inserted in different locations
		insertion_type: the location of insertion can be "assigned" by user or can be choose "random"
		insertion_locations: the list of locations to insert the AMR gene in case
			insertion_type = assigned
		amr_file:	the address of AMR file
		output_dir:	the path for the output directory
		metagenome_file: the name of the metagenome file
	Return:
		The list of name of genome files with AMR inserted in them and the address of metagenome_file
	"""
	#some validations
	if len(params.number_of_insertions) != len(params.ref_genome_files):
		logging.error("ERROR: Please specify the number of insertions for each reference genome.")
		#print("ERROR: Please specify the number of insertions for each reference genome.")
		sys.exit()
	if params.insertion_type == Insertion_type.assigned and \
		len(params.insertion_locations)!= sum(params.number_of_insertions):
		logging.error("ERROR: Please specify the location of insertion for all insertions OR \
			if you prefer them to be chosen randomely choose 1 for --insertion_type")
		#print("ERROR: Please specify the location of insertion for all insertions OR \
		#	if you prefer them to be chosen randomely choose 1 for --insertion_type")
		sys.exit()

	if not os.path.exists(output_dir):
		os.makedirs(output_dir)
	if (os.path.isfile(output_dir+metagenome_file)):
		logging.error('ERROR: A file named ' + metagenome_file + ' already exits in '+output_dir)
		#print('ERROR: A file named ' + metagenome_file + ' already exits in '+output_dir)
		sys.exit()

	insert_counter = 0
	genome_files = []
	amr_seq = retrieve_AMR(amr_file)
	for genome, insert_number in zip(ref_genome_files, number_of_insertions):
		if insertion_type.name == 'random':
			#find the number of lines in genome file
			lc = sum(1 for l in open(genome))
		for insertions in range(insert_number):
			if insertion_type.name == 'random':
				insert_location = random.randint(1,lc)
				new_file = insert_amr_in_location(genome, amr_seq, insert_location, output_dir)
			else:
				new_file = insert_amr_in_location(genome, amr_seq, insertion_locations[insert_counter], output_dir)
			genome_files.append(new_file)
			command = 'cat '+new_file+' >> ' + output_dir + metagenome_file
			os.system(command)
			insert_counter+=1
	return genome_files, output_dir + metagenome_file

def create_metagenome_file(ref_genome_files, output_dir, metagenome_file = 'metagenome.fasta'):
	"""
	To create a metagenome sample from some genomes
	Parameters:
		ref_genome_files:	the list of reference genomes in which AMR gene is supposed to be inserted
		output_dir:	the path for the output directory
		metagenome_file: the name of the metagenome file
	Return:
		The address of metagenome_file
	"""
	for genome in ref_genome_files:
		command = 'cat '+genome+' >> ' + output_dir + metagenome_file
		os.system(command)
	return output_dir + metagenome_file

def simulate_reads(metagenome_file, read_length, art_path, fcov = 20, sdev = 10):
	"""
	To use ART and simulate paired end reads for the metagenome sample
	Parameters:
		metagenome_file:	the address of the file containing the metagenome sample
		read_length:		the length of reads
		art_path:	the address of ART to run
		fcov: the fold of read coverage to be simulated or number of reads/read
			pairs generated for each amplicon
		sdev:	the standard deviation of DNA/RNA fragment size for paired-end
	Return:
		the address of paired end reads generated (read1 and read2)
	"""
	metagenome_file_name = os.path.splitext(metagenome_file)[0]
	if read_length == 150:
		command = art_path+' -ss HS25 -sam -i ' + metagenome_file + ' -p -l 150 \
			-f '+str(fcov)+' -m 500 -s '+str(sdev)+' -o ' + metagenome_file_name +'_'
	elif read_length ==250:
		command = art_path+' -ss MSv3 -sam -i ' + metagenome_file + ' -p -l 250 \
			-f '+str(fcov)+' -m 1000 -s '+str(sdev)+' -o ' + metagenome_file_name +'_'
	else:
		logging.error("ERROR: the read_length is not valid! In the current implementation it can\
		 be either 150 or 250!")
		# print("ERROR: the read_length is not valid! In the current implementation it can\
		#  be either 150 or 250!")
		sys.exit()
	logging.info("Running ART: "+command)
	os.system(command)

	read1 = metagenome_file_name + '_1.fq'
	read2 = metagenome_file_name + '_2.fq'
	return [read1, read2]

#def do_assembly(read1, read2, spades_path, output_dir, thread_num = 16, error_correction = True):
def do_assembly(reads, spades_path, output_dir, thread_num = 16, error_correction = True):
	"""
	To call MetaSPAdes and run assembly
	Parameters:
		reads:	the address of paired end reads either separated or interleaved
		spades_path:	the address of SPAdes to run from
		output_dir:		the path for the output directory
		thread_num:		the number of threads used in MetaSPAdes
		error_correction:to choose if Spades's error correction feature is used or not
	Return:
		the address of the file containing assembly graph (gfa_file) and
			contigs (contig_file)
	"""
	command = ""
	if len(reads)==2:
		command = 'python ' + spades_path + ' -1 '+reads[0]+' -2 '+reads[1]+' --meta '
	else:
		command = 'python ' + spades_path + ' --12 '+reads + ' --meta '
	if not error_correction:
		command+='--only-assembler '
	command += '--threads ' + str(thread_num) + ' -o '+ output_dir
	# if error_correction:
	# 	command = 'python ' + spades_path + ' -1 '+read1+' -2 '+read2+' --meta --threads ' \
	# 				+ str(thread_num) + ' -o '+ output_dir
	# else:
	# 	command = 'python ' + spades_path + ' -1 '+read1+' -2 '+read2+' --meta \
	# 	--only-assembler --threads ' + str(thread_num) + ' -o '+ output_dir
	logging.info("Running MetaSPAdes: "+command)
	os.system(command)

	#remove paths from GFA file
	gfa_file = output_dir + '/'+ ASSEMBLY_FILE
	command = "sed -i '/^P/d' " + gfa_file
	os.system(command)

	contig_file = output_dir + '/' + CONTIG_FILE

	return gfa_file, contig_file

def run_RGI(input_file, output_dir, seq_description, include_loose = False):
	"""
	To run RGI and annotate AMRs in the sequence
	# To ensure consistency between Prokka and RGI findings, we annotate found proteins
	# by Prokka (instead of annotationg DNA sequences from scratch)
	Parameters:
		input_file: the file contating proteins annotated by Prokka
		output_dir:  the path for the output directory
		seq_description: a small description of the sequence used for naming
		include_loose: Whether to include loose annotations
	Return:
		the list of extracted annotation information for the sequence
	"""
	rgi_dir = output_dir +"rgi_dir"
	if not os.path.exists(rgi_dir):
		os.makedirs(rgi_dir)
	output_file_name = rgi_dir +"/rgi_output_"+seq_description+"_"+datetime.datetime.now().strftime('%Y-%m-%d_%H-%M')
	#remove any potential * from the sequence
	command = "sed -i 's/*//g' " + input_file
	os.system(command)
	command = "rgi main --input_sequence " + input_file + " --output_file " +\
		output_file_name +" --input_type protein --clean --exclude_nudge"
	if include_loose:
		command+=" --include_loose"
	#rgi main --input_sequence  <seq_file_name> --output_file <output_file_name> --clean --include_loose --exclude_nudge --low_quality
	os.system(command)
	seq_info_list = []
	if os.path.isfile(output_file_name + '.txt'):
		with open(output_file_name + '.txt', newline = '') as rgi_file:
			rgi_reader = csv.reader(rgi_file, delimiter='\t')
			next(rgi_reader)
			for row in rgi_reader:
				seq_info = {'ORF_ID':row[0], 'gene':row[8].strip(),
				'prediction_type':row[5].strip(), 'best_identities':float(row[9]),
				'family':row[16].strip()}
				seq_info_list.append(seq_info)
	else:
		logging.error("ERROR: RGI didn't run successfully!")
		#print("ERROR: RGI didn't run successfully!")
		sys.exit()

	return seq_info_list


def annotate_sequence(seq, seq_description, output_dir, prokka_prefix, use_RGI = True,\
						RGI_include_loose = False):
	"""
	To run Prokka for a sequence and extract required information from its
		generated output files
	Parameters:
		seq:	the sequence to be annotated
		seq_description: a small description of the sequence used for naming
		output_dir:  the path for the output directory
		use_RGI:	RGI annotations incorporated for AMR annotation
		prokka_prefix: to run prokka via docker or conda or any other source properly
	Return:
		the list of extracted annotation information for the sequence
	"""
	#write the sequence into a temporary file
	seq_file_name = create_fasta_file(seq, output_dir, file_name='temp_'+seq_description)

	prokka_dir = output_dir+'prokka_dir_'+seq_description+'_'+datetime.datetime.now().strftime('%Y-%m-%d_%H-%M')
	prefix_name = 'mygenome_'+seq_description
	# if os.path.exists(prokka_dir):
	# 	try:
	# 		shutil.rmtree(prokka_dir)
	# 	except OSError as e:
	# 		logging.error("Error: %s - %s." % (e.filename, e.strerror))
			#print("Error: %s - %s." % (e.filename, e.strerror))
	command = prokka_prefix + 'prokka --metagenome --outdir '+\
		prokka_dir+' --prefix '+ prefix_name+' --fast --notrna '+seq_file_name
	os.system(command)

	RGI_output_list = None
	if use_RGI:
		RGI_output_list = run_RGI(prokka_dir+'/'+prefix_name+'.faa', output_dir, seq_description,
								RGI_include_loose)

	#Go over Prokka's output files and extract required information
	seq_info = []
	with open(prokka_dir+'/'+prefix_name+'.tsv') as tsvfile:
		reader = csv.reader(tsvfile, delimiter='\t')
		#skip the header
		next(reader)
		for row in reader:
			mygene = row[3].strip()
			split_gene = mygene.split('_')
			if len(split_gene)==2 and split_gene[1].isnumeric():
				mygene = split_gene[0]
			gene_info = {'locus_tag':row[0].strip(), 'gene':mygene,'length':row[2].strip(),
						'product':row[6].strip(),'start_pos':None, 'end_pos':None,
						'prokka_gene_name':mygene, 'RGI_prediction_type':None,
						'coverage':None, 'family': None, 'seq_value': seq[:-1],
						'seq_name':None, 'target_amr': None}
			seq_info.append(gene_info)
	counter = 0
	with open(prokka_dir+'/'+prefix_name+'.tbl', 'r') as read_obj:
		for line in read_obj:
			if line[0].isdigit():
				cells = line.split('\t')
				seq_info[counter]['start_pos'] = int(cells[0])
				seq_info[counter]['end_pos'] = int(cells[1])
				counter+=1

	#incorporate RGI findings into Prokka's
	if RGI_output_list:
		for item in RGI_output_list:
			for gene_info in seq_info:
				if item['ORF_ID'].split(' ')[0]==gene_info['locus_tag']:
					gene_info['gene'] = item['gene']
					gene_info['RGI_prediction_type'] = item['prediction_type']
					gene_info['family'] = item['family']
					break

	if os.path.isfile(seq_file_name):
		os.remove(seq_file_name)

	# #remove temporary files and folder
	# os.remove(seq_file_name)
	# try:
	# 	shutil.rmtree(prokka_dir)
	# except OSError as e:
	# 	print("Error: %s - %s." % (e.filename, e.strerror))
	return seq_info

def unnamed_genes_are_siginificantly_similar(gene_info1, gene_info2, threshold = 90):
	"""
	"""
	if gene_info1['gene']!='' or gene_info2['gene']!='':
		return False
	start1, end1 = min(gene_info1['start_pos'], gene_info1['end_pos']), max(gene_info1['start_pos'], gene_info1['end_pos'])
	seq1 = gene_info1['seq_value'][start1-1:end1-1]
	start2, end2 = min(gene_info2['start_pos'], gene_info2['end_pos']), max(gene_info2['start_pos'], gene_info2['end_pos'])
	seq2 = gene_info2['seq_value'][start2-1:end2-1]
	return compare_two_sequences(seq1, seq2, '', threshold)

def seqs_annotation_are_identical(seq_info1, seq_info2, threshold = 90):
	"""
	"""
	if len(seq_info1)==len(seq_info2):
		identical_rows = 0
		for i, gene_info1 in enumerate(seq_info1):
			gene_info2 = seq_info2[i]
			if (gene_info1['gene']==gene_info2['gene'] and gene_info1['gene']!='') or\
				(gene_info1['gene']==gene_info2['gene'] and\
				unnamed_genes_are_siginificantly_similar(gene_info1, gene_info2, threshold) ):
				identical_rows+=1
		if identical_rows == len(seq_info1):
			return True
	return False

def seq_annotation_already_exist(seq_info_list, all_seq_info_lists):
	"""
	To check if annotations found for the new sequence have already exists in the
	list of annotations extracted from other sequences.
	Parameters:
		seq_info_list: list of annotations of the new sequence
		all_seq_info_lists: list of annotations of the other sequences annotated so far
	Return:
		True if seq_info_list is identical to another list in all_seq_info_lists in
		terms of 'gene', 'length','start_pos' and 'end_pos' for all thei items
	"""
	found = False
	for seq_list in all_seq_info_lists:
		if seqs_annotation_are_identical(seq_info_list, seq_list, 100):
		# if len(seq_list)==len(seq_info_list):
		# 	identical_rows = 0
		# 	for i, seq_info in enumerate(seq_info_list):
		# 		if (seq_info['gene']==seq_list[i]['gene'] or\
		# 			(seq_info['gene']=='' and (seq_list[i]['gene']=='UNKNOWN' or\
		# 										seq_list[i]['gene']==''))) and\
		# 			seq_info['length']==seq_list[i]['length'] and\
		# 		 	seq_info['start_pos']==seq_list[i]['start_pos'] and\
		# 			seq_info['end_pos']==seq_list[i]['end_pos']:
		# 			identical_rows+=1
		# 	if identical_rows == len(seq_list):
			found = True
			break

	return found

def similar_seq_annotation_already_exist(seq_info_list, all_seq_info_lists, threshold = 90):
	"""
	To check if annotations found for the new sequence have already exists in the
	list of annotations extracted from other sequences.
	Parameters:
		seq_info_list: list of annotations of the new sequence
		all_seq_info_lists: list of annotations of the other sequences annotated so far
	Return:
		True if seq_info_list is identical to another list in all_seq_info_lists in
		terms of 'gene', 'length','start_pos' and 'end_pos' for all thei items
	"""
	found = False
	for seq_list in all_seq_info_lists:
		if seqs_annotation_are_identical(seq_info_list, seq_list, threshold):
		# if len(seq_list)==len(seq_info_list):
		# 	similar_rows = 0
		# 	for i, seq_info in enumerate(seq_info_list):
		# 		if (seq_info['gene']==seq_list[i]['gene'] and seq_info['gene']!='') or\
		# 			(seq_info['gene']==seq_list[i]['gene'] and\
		# 			unnamed_genes_are_siginificantly_similar(seq_info, seq_list[i], threshold = 90)):
		# 			#(seq_info['gene']==seq_list[i]['gene'] and seq_info['length']==seq_list[i]['length']):
		# 		# if (seq_info['gene']==seq_list[i]['gene'] or\
		# 		# 	(seq_info['gene']=='' and (seq_list[i]['gene']=='UNKNOWN' or\
		# 		# 								seq_list[i]['gene']==''))):
		# 		# 	if seq_info['length']==seq_list[i]['length'] or\
		# 		#  	(seq_info['start_pos']==seq_list[i]['start_pos'] and\
		# 		# 	seq_info['end_pos']==seq_list[i]['end_pos']):
		# 			similar_rows+=1
		# 	if similar_rows == len(seq_list):
			found = True
			break

	return found


def read_path_info_file(path_info_file):
	"""
	"""
	seq_info_list = []
	seq_info = []
	with open(path_info_file, 'r') as myfile:
		myreader = DictReader(myfile)
		old_seq = ''
		for row in myreader:
			node_info = {'node':row['node'], 'coverage':float(row['coverage']),
			'start':int(row['start']), 'end':int(row['end'])}
			cur_seq = row['sequence']
			if cur_seq!=old_seq:
				if (seq_info):
					seq_info_list.append(seq_info)
				seq_info = []
				old_seq = cur_seq
			seq_info.append(node_info)
		seq_info_list.append(seq_info)
	return seq_info_list

def find_gene_coverage(seq_info_list, path_info):
	"""
	"""
	coverage_list = []
	for seq_info in seq_info_list:
		sum_coverage = 0
		#minus 1 because I guess in what prokka returns the sequence starts from position 1
		start, end = seq_info['start_pos']-1, seq_info['end_pos']-1
		if start > end:
			start, end = end, start
		found = False
		for j, node_info in enumerate(path_info):
			if node_info['start']<= start and node_info['end']>= start:
				found = True
				#the annotated gene is represented by a single node
				if node_info['end'] >= end:
					sum_coverage = (end - start + 1) * node_info['coverage']
				else:
					#calculate the length of node representing the gene
					sum_coverage = (node_info['end'] - start + 1)*node_info['coverage']
					n_index = j + 1
					assert n_index < len(path_info), "wrong index calculated for path_info!!"
					while (n_index < len(path_info) and path_info[n_index]['start']<=end):
						if path_info[n_index]['end'] >= end:
							sum_coverage+=(end-path_info[n_index]['start']+1)*path_info[n_index]['coverage']
							break
						else:
							sum_coverage+=(path_info[n_index]['end']-path_info[n_index]['start']+1)*path_info[n_index]['coverage']
						n_index+=1
				coverage_list.append(sum_coverage / (end-start+1))
				break
		if not found:
			print("ERROR: no nodes were found for this gene!!!")
			import pdb; pdb.set_trace()
			sys.exit()
	return coverage_list

def this_is_target_amr_row(amr_name, gene, heads, member_lists):
	"""
	"""
	gene = gene.strip().replace(' ','_')
	gene_name = ''.join(e for e in gene if e.isalpha() or e.isnumeric() or e=='_' or e=='-')
	if amr_name.strip()==gene_name.strip():
		return True
	#else:
	#found = False
	for i, item in enumerate(heads):
		trimmed_item = ''.join(e for e in item if e.isalpha() or e.isnumeric() or e=='_' or e=='-')
		if trimmed_item==amr_name:
			#found = True
			for member in member_lists[i]:
				if member.strip()==gene.strip():
					return True
	# if not found:
	# 	logging.info("ERROR: not able to find "+amr_name+" in overlap file as a head!!!")
	# 	sys.exit()
	return False

def find_target_amr_in_seqvalue_and_return_coverage(seq_info):
	"""
	"""
	error = True
	amr_coverage = 0
	#find the indeces of lower case string (amr sequence) in the extracted sequence
	sequence = seq_info[0]['seq_value']
	amr_start = -1
	amr_end = -1
	index = 0
	while index < len(sequence):
		if sequence[index].islower() and amr_start==-1:
			amr_start = index
		elif sequence[index].isupper() and amr_start>-1:
			amr_end=index-1
			break
		index+=1
	#import pdb; pdb.set_trace()
	#find the gene has the most overlap with the found range
	#we look at all found genes that their overlap with the sequence is more than initial value of most_overlap and chose the one with max
	most_overlap = 50
	amr_index = -1
	for i, gene_info in enumerate(seq_info):
		start, end = min(gene_info['start_pos'], gene_info['end_pos']), max(gene_info['start_pos'], gene_info['end_pos'])
		if end<amr_start or start> amr_end:
			continue
		else:
			# added by 1 because in string indecesstarts from 0
			diff = max((amr_start+1-start), 0)+max((end - (amr_end+1)), 0)
			if ((1-(float(diff)/(end-start)))*100)>most_overlap:
				most_overlap = (1-(float(diff)/(end-start)))*100
				amr_coverage = gene_info['coverage']
				amr_index = i
				error = False
	return amr_coverage, amr_index, error

def check_coverage_consistency_remove_rest_seq(seq_info_list_input, ref_seq_info_list,
								output_dir, coverage_thr, amr_name, annotate_dir):
	"""
	"""
	seq_info_list = copy.deepcopy(seq_info_list_input)
	#extract the list of AMR groups
	overlap_file_name = output_dir+AMR_DIR_NAME+AMR_OVERLAP_FILE
	heads = []
	member_lists = []
	with open(overlap_file_name, 'r') as read_obj:
		for line in read_obj:
			if line[:-1].endswith(':'):
				heads.append(line[:-2])
			else:
				members = line[:-1].split(', ')
				member_lists.append(members)
	#extract amr info
	amr_coverages = []
	amr_indeces = []
	for seq_info in seq_info_list:
		found_amr = False
		for gene_counter, gene_info in enumerate(seq_info):
			if gene_info['coverage'] is None:
				logging.info("Coverage information are not available for "+amr_name)
				return "", 0
			coverage = round(gene_info['coverage'], 2)
			#find the target AMR gene
			#target_amr =  this_is_target_amr_row(amr_name, gene_info['gene'], heads, member_lists)
			#if target_amr:
			if gene_info['target_amr']=='yes':
				amr_coverages.append(coverage)
				amr_indeces.append(gene_counter)
				found_amr = True
				break
		if not found_amr:
			amr_coverage, amr_index, error = find_target_amr_in_seqvalue_and_return_coverage(seq_info)
			if error:
				logging.info("ERROR: no target amr was found for "+ str(seq_info)+" regarding "+amr_name)
				sys.exit()
			else:
				amr_coverages.append(amr_coverage)
				amr_indeces.append(amr_index)
	if len(amr_coverages)!=len(seq_info_list):
		logging.error("ERROR: inconsistency between the number of sequences and found amr-coverages!")
		import pdb; pdb.set_trace()
	#check coverage consistency by comparing its coverage with AMR coverage
	# and remove genes with inconsistent coverage and whatever comes before them if upstream OR after them if downstream
	remained_seqs = []
	for i, seq_info in enumerate(seq_info_list):
		#find the genes need to be removed
		to_be_removed_genes=[]
		for j, gene_info in enumerate(seq_info):
			if abs(gene_info['coverage'] - amr_coverages[i])>coverage_thr:
				if j<amr_indeces[i]:
					for k in range(j+1):
						if k not in to_be_removed_genes:
							to_be_removed_genes.append(k)
				elif j>amr_indeces[i]:
					for k in range(j, len(seq_info)):
						if k not in to_be_removed_genes:
							to_be_removed_genes.append(k)
					break
		for j in reversed(range(len(seq_info))):
			if j in to_be_removed_genes:
				del seq_info[j]
		#check if the remained sequence already exists in the seq_info_list
		if seq_info and not similar_seq_annotation_already_exist(seq_info, remained_seqs):
			remained_seqs.append(seq_info)

	#Initialize coverage file
	coverage_annotation = annotate_dir+'coverage_annotation_'+str(coverage_thr)+'_'+amr_name+'.csv'
	with open(coverage_annotation,'w') as fd:
		writer = csv.writer(fd)
		writer.writerow(['seq_name', 'seq_value', 'seq_length', 'gene',
							'coverage', 'length', 'start_pos', 'end_pos', 'target_amr'])
		#copy reference genomes info from annotation_info file
		for ref_info in ref_seq_info_list:
			for row in ref_info:
				writer.writerow([row['seq_name'], row['seq_value'], len(row['seq_value']),
								row['gene'], '', row['length'], row['start_pos'],
								row['end_pos'], row['target_amr']])
		#write extracted sequences with consistent coverage
		for seq_info in remained_seqs:
			for gene_info in seq_info:
				writer.writerow([gene_info['seq_name'], gene_info['seq_value'],
							len(gene_info['seq_value']),
							gene_info['gene'], gene_info['coverage'],
							gene_info['length'], gene_info['start_pos'],
							gene_info['end_pos'], gene_info['target_amr']])

	# #check its coverage and compare with other genes' coverage in the sequence
	# to_be_removed_seqs = []
	# for i, seq_info in enumerate(seq_info_list):
	# 	for gene_info in seq_info:
	# 		if abs(gene_info['coverage'] - amr_coverages[i])>coverage_thr:
	# 			#if only one gene is not consistent with AMR remove the entire sequence
	# 			to_be_removed_seqs.append(i)
	# 			break
	return coverage_annotation, len(remained_seqs)

def check_coverage_consistency(annotation_info, output_dir, coverage_thr, amr_name,
								annotate_dir):
	"""
	"""
	#extract the list of AMR groups
	overlap_file_name = output_dir+AMR_DIR_NAME+AMR_OVERLAP_FILE
	heads = []
	member_lists = []
	with open(overlap_file_name, 'r') as read_obj:
		for line in read_obj:
			if line[:-1].endswith(':'):
				heads.append(line[:-2])
			else:
				members = line[:-1].split(', ')
				member_lists.append(members)
	#Read the csv file and store extraction info
	seq_info_list = []
	seq_info = []
	amr_coverages = []
	amr_coverage = -1
	with open(annotation_info, 'r') as myfile:
		myreader = DictReader(myfile)
		old_seq = ''
		for row in myreader:
			if row['seq_name'].startswith('extracted'):
				if row['coverage']=='':
					logging.info("Coverage information are not available for "+amr_name)
					return "", 0
				#else:
				coverage = round(float(row['coverage']), 2)
				#find the target AMR gene
				if this_is_target_amr_row(amr_name, row['gene'], heads, member_lists):
					amr_coverage = coverage
				gene_info = {'seq_name':row['seq_name'], 'seq_value':row['seq_value'],
						'seq_length':row['seq_length'],
						'gene':row['gene'], 'coverage':coverage, 'length':row['length'],
						'start_pos':int(row['start_pos']), 'end_pos':int(row['end_pos'])}
				cur_seq = row['seq_name']
				if cur_seq!=old_seq:
					if (seq_info):
						seq_info_list.append(seq_info)
						if amr_coverage!=-1:
							amr_coverages.append(amr_coverage)
						else:
							amr_coverage, _, error = find_target_amr_in_seqvalue_and_return_coverage(seq_info)
							if error:
								logging.info("ERROR: no target amr was found for "+ str(seq_info)+" regarding "+amr_name)
								sys.exit()
							else:
								amr_coverages.append(amr_coverage)
						amr_coverage = -1
						seq_info = []
					old_seq = cur_seq
				seq_info.append(gene_info)
		seq_info_list.append(seq_info)
		if amr_coverage!=-1:
			amr_coverages.append(amr_coverage)
		else:
			amr_coverage, _, error = find_target_amr_in_seqvalue_and_return_coverage(seq_info)
			if error:
				logging.info("ERROR: no target amr was found for "+ str(seq_info)+" regarding "+amr_name)
				sys.exit()
			else:
				amr_coverages.append(amr_coverage)

	if len(amr_coverages)!=len(seq_info_list):
		logging.error("ERROR: inconsistency between the number of sequences and found amr-coverages!")
		import pdb; pdb.set_trace()
	#check coverage consistency by comparing its coverage with AMR coverage
	# and remove genes with inconsistent coverage
	remained_seqs = []
	for i, seq_info in enumerate(seq_info_list):
		#find the genes need to be removed
		to_be_removed_genes=[]
		for j, gene_info in enumerate(seq_info):
			if abs(gene_info['coverage'] - amr_coverages[i])>coverage_thr:
				to_be_removed_genes.append(j)
		for j in reversed(range(len(seq_info))):
			if j in to_be_removed_genes:
				del seq_info[j]
		#check if the remained sequence already exists in the seq_info_list
		if not similar_seq_annotation_already_exist(seq_info, remained_seqs):
			remained_seqs.append(seq_info)

	#Initialize coverage file
	coverage_annotation = annotate_dir+'coverage_annotation.csv'
	with open(coverage_annotation,'w') as fd:
		writer = csv.writer(fd)
		writer.writerow(['seq_name', 'seq_value', 'seq_length', 'gene',
							'coverage', 'length', 'start_pos', 'end_pos'])
		#copy reference genomes info from annotation_info file
		with open(annotation_info, 'r') as myfile:
			myreader = DictReader(myfile)
			old_seq = ''
			for row in myreader:
				if not row['seq_name'].startswith('extracted'):
					writer.writerow([row['seq_name'], row['seq_value'], row['seq_length'],
									row['gene'],
									'', row['length'], row['start_pos'], row['end_pos']])
		#write extracted sequences with consistent coverage
		for seq_info in remained_seqs:
			for gene_info in seq_info:
				writer.writerow([gene_info['seq_name'], gene_info['seq_value'],
							len(gene_info['seq_value']),
							gene_info['gene'], gene_info['coverage'],
							gene_info['length'], gene_info['start_pos'],
							gene_info['end_pos']])

	# #check its coverage and compare with other genes' coverage in the sequence
	# to_be_removed_seqs = []
	# for i, seq_info in enumerate(seq_info_list):
	# 	for gene_info in seq_info:
	# 		if abs(gene_info['coverage'] - amr_coverages[i])>coverage_thr:
	# 			#if only one gene is not consistent with AMR remove the entire sequence
	# 			to_be_removed_seqs.append(i)
	# 			break
	return coverage_annotation, len(remained_seqs), amr_coverages

def split_up_down_info(sequence, seq_info):
	"""
	"""
	amr_start = -1
	amr_end = -1
	index = 0
	up_info = []
	down_info = []
	while index < len(sequence):
		if sequence[index].islower() and amr_start==-1:
			amr_start = index
		elif sequence[index].isupper() and amr_start>-1:
			amr_end=index-1
			break
		index+=1
	#import pdb;pdb.set_trace()
	#find the gene has the most overlap with the found range
	overlap_thr = 50
	found = False
	amr_info = []
	for gene_info in seq_info:
		start, end = min(gene_info['start_pos'], gene_info['end_pos']), max(gene_info['start_pos'], gene_info['end_pos'])
		if end<amr_start:
			up_info.append(gene_info)
		elif start> amr_end:
			down_info.append(gene_info)
		else:
			# added by 1 because in string indecesstarts from 0
			diff = max((amr_start+1-start), 0)+max((end - (amr_end+1)), 0)
			if ((1-(float(diff)/(end-start)))*100)>overlap_thr:
				found = True
				gene_info['target_amr'] = 'yes'
				amr_info = gene_info
			elif start<amr_start:
				up_info.append(gene_info)
			else:
				down_info.append(gene_info)

	return found, amr_info, up_info, down_info, seq_info

def exists_in_gene_list(gene_info_list, gene_info):
	"""
	"""
	for item in gene_info_list:
		if gene_info['gene']!='' and item['gene']==gene_info['gene']:
			return True
		elif gene_info['gene']=='' and item['gene']==gene_info['gene']:
			return unnamed_genes_are_siginificantly_similar(item, gene_info, threshold = 90)
			# if gene_info['start_pos']==item['start_pos'] and gene_info['end_pos'] == item['end_pos']:
			# 	return True
			# #for now we assume even if their length is equal they are identical!!!
			# elif gene_info['length'] == item['length']:
			# 	return True
	return False

def extract_up_down_from_csv_file(seq_info):
	"""
	"""
	up_info = []
	down_info = []
	amr_info = []
	amr_found = False
	for gene_info in seq_info:
		if gene_info['target_amr']=='yes':
			amr_info = gene_info
			amr_found = True
		elif not amr_found:
			up_info.append(gene_info)
		else:
			down_info.append(gene_info)

	return amr_found, up_info, down_info, amr_info

def evaluate_sequences_up_down_based_on_coverage(amr_name, coverage_annotation, summary_file,
													ref_up_info_list, ref_down_info_list, ref_amr_info_list):
	"""
	"""
	#Read  coverage_annotation csv file and store info
	seq_info = []
	#seq_info_list = []
	up_info_list = []
	down_info_list = []
	amr_info_list = []
	sequence = ''
	with open(coverage_annotation, 'r') as myfile:
		myreader = DictReader(myfile)
		old_seq = ''
		for row in myreader:
			if row['seq_name'].startswith('extracted'):
				if row['coverage']=='':
					logging.error("ERROR: Coverage information are not available for "+amr_name)
					sys.exit()
				#else:
				gene_info = {'seq_name':row['seq_name'], 'seq_value':row['seq_value'],
				 			'gene':row['gene'], 'length':row['length'],
							'start_pos':int(row['start_pos']),'end_pos':int(row['end_pos']),
							'target_amr':row['target_amr']}
				cur_seq = row['seq_name']
				if cur_seq!=old_seq:
					if (seq_info):
						#amr_found, amr_info, up_info, down_info, seq_info = split_up_down_info(sequence, seq_info)
						amr_found, up_info, down_info, amr_info = extract_up_down_from_csv_file(seq_info)
						if amr_found:
							amr_info_list.append(amr_info)
							if up_info and not similar_seq_annotation_already_exist(up_info, up_info_list):
								up_info_list.append(up_info)
							if down_info and not similar_seq_annotation_already_exist(down_info, down_info_list):
								down_info_list.append(down_info)
					seq_info = []
					old_seq = cur_seq
				seq_info.append(gene_info)
				sequence = row['seq_value']
		#seq_info_list.append(seq_info)
		#amr_found, amr_info, up_info, down_info, seq_info = split_up_down_info(sequence, seq_info)
		amr_found, up_info, down_info, amr_info = extract_up_down_from_csv_file(seq_info)
		if amr_found:
			amr_info_list.append(amr_info)
			if up_info and not similar_seq_annotation_already_exist(up_info, up_info_list):
				up_info_list.append(up_info)
			if down_info and not similar_seq_annotation_already_exist(down_info, down_info_list):
				down_info_list.append(down_info)

	ref_len =  len(ref_up_info_list)+len(ref_down_info_list)
	if ref_len==0 and len(ref_amr_info_list)==0:
		with open(summary_file,'a') as fd:
			writer = csv.writer(fd)
			writer.writerow([amr_name, 0, 0, 0, len(up_info_list)+len(down_info_list),-1, -1])
		return -1, -1
	#find the number of unique true-positives, all false positives, total found cases, all unique true cases
	unique_tp = 0
	for ref_info in ref_up_info_list:
		for seq_info in up_info_list:
			if seqs_annotation_are_identical(ref_info, seq_info):
				unique_tp+=1
				break
	for ref_info in ref_down_info_list:
		for seq_info in down_info_list:
			if seqs_annotation_are_identical(ref_info, seq_info):
				unique_tp+=1
				break
	sensitivity = 1 if ref_len==0 else round(float(unique_tp)/ref_len, 2)

	#We probably don't need to calculate false positive for this case because annotations are unique anyway
	false_positive = 0
	found_cases_len = 0
	for seq_info in up_info_list:
		if seq_info:
			found_cases_len+=1
			found_identical_annotation = False
			for ref_info in ref_up_info_list:
				if seqs_annotation_are_identical(ref_info, seq_info):
					found_identical_annotation = True
					break
			if not found_identical_annotation:
				false_positive+=1
	for seq_info in down_info_list:
		if seq_info:
			found_cases_len+=1
			found_identical_annotation = False
			for ref_info in ref_down_info_list:
				if seqs_annotation_are_identical(ref_info, seq_info):
					found_identical_annotation = True
					break
			if not found_identical_annotation:
				false_positive+=1
	#found_cases_len = len(up_info_list)+len(down_info_list)
	if found_cases_len == 0 and ref_len == 0:
		precision = 1
	else:
		precision = 0 if found_cases_len==0 else round(1 - float(false_positive)/found_cases_len, 2)

	with open(summary_file,'a') as fd:
		writer = csv.writer(fd)
		writer.writerow([amr_name, unique_tp, false_positive, ref_len, found_cases_len,
						sensitivity, precision])

	return sensitivity, precision

def evaluate_sequences_gene_based_on_coverage(amr_name, coverage_annotation, summary_file,
											annotation_file, amr_coverages):
	"""
	"""
	#Read  coverage_annotation csv file and store info
	ref_info_list = []
	seq_info = []
	seq_info_list = []
	with open(coverage_annotation, 'r') as myfile:
		myreader = DictReader(myfile)
		old_seq = ''
		for row in myreader:
			if not row['seq_name'].startswith('extracted'):
				ref_info = {'seq_name':row['seq_name'], 'seq_value':row['seq_value'],
				 			'gene':row['gene'],'length':row['length'],
							'start_pos':int(row['start_pos']), 'end_pos':int(row['end_pos'])}
				#find the number of different genes in ref_genomes
				if not exists_in_gene_list(ref_info_list, ref_info):
					ref_info_list.append(ref_info)
			else:
				if row['coverage']=='':
					logging.error("ERROR: Coverage information are not available for "+amr_name)
					sys.exit()
				#else:
				gene_info = {'seq_name':row['seq_name'], 'gene':row['gene'],
						'length':row['length'],'start_pos':int(row['start_pos']),
						'end_pos':int(row['end_pos'])}
				cur_seq = row['seq_name']
				if cur_seq!=old_seq:
					if (seq_info):
						seq_info_list.append(seq_info)
					seq_info = []
					old_seq = cur_seq
				seq_info.append(gene_info)
		seq_info_list.append(seq_info)

	if not ref_info_list:
		return -1, -1
	#find the number of unique true-positives, all false positives, total found cases, all unique true cases
	true_positive = 0
	for ref_info in ref_info_list:
		found = False
		for seq_info in seq_info_list:
			found = exists_in_gene_list(seq_info, ref_info)
			if found:
				true_positive+=1
				break
	sensitivity = round(float(true_positive)/len(ref_info_list), 2)
	false_positive = 0
	total = 0
	for seq_info in seq_info_list:
		total+=len(seq_info)
		for gene_info in seq_info:
			if not exists_in_gene_list(ref_info_list, gene_info):
				false_positive+=1
	precision = round(1 - float(false_positive)/total, 2)

	#Read  annotation csv file and store info
	seq_info = []
	seq_info_list2 = []
	with open(annotation_file, 'r') as myfile:
		myreader = DictReader(myfile)
		old_seq = ''
		for row in myreader:
			if row['seq_name'].startswith('extracted'):
				coverage = round(float(row['coverage']), 2)
				gene_info = {'seq_name':row['seq_name'], 'gene':row['gene'],
						'coverage':coverage, 'length':row['length'],
						'start_pos':int(row['start_pos']), 'end_pos':int(row['end_pos'])}
				cur_seq = row['seq_name']
				if cur_seq!=old_seq:
					if (seq_info):
						seq_info_list2.append(seq_info)
						seq_info = []
					old_seq = cur_seq
				seq_info.append(gene_info)
		seq_info_list2.append(seq_info)
	if len(amr_coverages)!=len(seq_info_list2):
		print("ERROR: inconsistency between length of amr_coverages and seq_info_list2")
		import pdb; pdb.set_trace()
	#find the min and max coverage for false and true cases respectively
	min_false = 1000
	found_false = False
	max_true = 0
	for i, seq_info in enumerate(seq_info_list2):
		for gene_info in seq_info:
			if exists_in_gene_list(ref_info_list, gene_info):
				if abs(gene_info['coverage']-amr_coverages[i])> max_true:
					max_true = round(float(abs(gene_info['coverage']-amr_coverages[i])),2)
			elif abs(gene_info['coverage']-amr_coverages[i]) < min_false:
				found_false = True
				min_false = round(float(abs(gene_info['coverage']-amr_coverages[i])), 2)

	with open(summary_file,'a') as fd:
		writer = csv.writer(fd)
		writer.writerow([amr_name, true_positive, false_positive, len(ref_info_list), total,
						sensitivity, precision, max_true, min_false if found_false else ''])

	return sensitivity, precision

def extract_seq_annotation(annotate_dir, prokka_prefix, use_RGI, RGI_include_loose,
								seq_pair):
	"""
	"""
	counter, ext_seq = seq_pair
	seq_description = 'extracted'+str(counter)
	seq_info_list = annotate_sequence(ext_seq, seq_description, annotate_dir+'/',
									prokka_prefix, use_RGI, RGI_include_loose)
	return seq_info_list

def add_RGI_loose(sequence, seq_info, up_info_list, down_info_list):
	"""
	To correct the annotation of extracted sequences
	if the prokka annotation for a part of the extracted sequence and ref sequence
	are the same and while there is an RGI annotation for that part of ref sequence,
	there is no RGI annotation for extracted sequence, we add ref RGI annotation to
	the extracted sequence annotation and consider this as an RGI loose hit!
	"""
	amr_start = -1
	amr_end = -1
	index = 0
	#up_info = []
	#down_info = []
	while index < len(sequence):
		if sequence[index].islower() and amr_start==-1:
			amr_start = index
		elif sequence[index].isupper() and amr_start>-1:
			amr_end=index-1
			break
		index+=1
	#find the gene has the most overlap with the found range
	overlap_thr = 50
	found = False
	#amr_info = []
	for gene_info in seq_info:
		start, end = min(gene_info['start_pos'], gene_info['end_pos']), max(gene_info['start_pos'], gene_info['end_pos'])
		if end<amr_start:
			#there is no RGI prediction for this gene
			if gene_info['gene']==gene_info['prokka_gene_name'] and gene_info['gene']!='':
				for ref_info in up_info_list:
					for ref_gene_info in ref_info:
						if ref_gene_info['prokka_gene_name']==gene_info['prokka_gene_name'] and\
							ref_gene_info['prokka_gene_name']!=ref_gene_info['gene']:
							gene_info['gene']=ref_gene_info['gene']
							gene_info['RGI_prediction_type'] = 'loose'
							break
					else:
						continue
					break
			#up_info.append(gene_info)
		elif start> amr_end:
			#there is no RGI prediction for this gene
			if gene_info['gene']==gene_info['prokka_gene_name'] and gene_info['gene']!='':
				for ref_info in down_info_list:
					for ref_gene_info in ref_info:
						if ref_gene_info['prokka_gene_name']==gene_info['prokka_gene_name'] and\
							ref_gene_info['prokka_gene_name']!=ref_gene_info['gene']:
							gene_info['gene']=ref_gene_info['gene']
							gene_info['RGI_prediction_type'] = 'loose'
							break
					else:
						continue
					break
			#down_info.append(gene_info)
		else:
			# added by 1 because in string indecesstarts from 0
			diff = max((amr_start+1-start), 0)+max((end - (amr_end+1)), 0)
			if ((1-(float(diff)/(end-start)))*100)>overlap_thr:
				found = True
				#amr_info = gene_info
			elif start<amr_start:
				#there is no RGI prediction for this gene
				if gene_info['gene']==gene_info['prokka_gene_name'] and gene_info['gene']!='':
					for ref_info in up_info_list:
						for ref_gene_info in ref_info:
							if ref_gene_info['prokka_gene_name']==gene_info['prokka_gene_name'] and\
								ref_gene_info['prokka_gene_name']!=ref_gene_info['gene']:
								gene_info['gene']=ref_gene_info['gene']
								gene_info['RGI_prediction_type'] = 'loose'
								break
						else:
							continue
						break
				#up_info.append(gene_info)
			else:
				#there is no RGI prediction for this gene
				if gene_info['gene']==gene_info['prokka_gene_name'] and gene_info['gene']!='':
					for ref_info in down_info_list:
						for ref_gene_info in ref_info:
							if ref_gene_info['prokka_gene_name']==gene_info['prokka_gene_name'] and\
								ref_gene_info['prokka_gene_name']!=ref_gene_info['gene']:
								gene_info['gene']=ref_gene_info['gene']
								gene_info['RGI_prediction_type'] = 'loose'
								break
						else:
							continue
						break
				#down_info.append(gene_info)
	return seq_info
	#return found, amr_info, up_info, down_info, seq_info

def add_coverage_to_info(seq_info, up_info, down_info, amr_info, coverage_list):
	"""
	"""
	for j, gene_info in enumerate(seq_info):
		coverage = coverage_list[j]
		gene_info['coverage'] = coverage
		if j<len(up_info):
			up_info[j]['coverage'] = coverage
		elif j<len(up_info)+1:
			amr_info['coverage'] = coverage
		else:
			down_info[j-len(up_info)-1]['coverage'] = coverage

	return seq_info, up_info, down_info, amr_info

def test_function():
	"""
	"""
	#import pdb; pdb.set_trace()
	amr_seq = retrieve_AMR('411.fasta')
	annotate_dir = 'Experiments/1_1_1_limit/test_411/'
	os.makedirs(annotate_dir)
	annotation_detail_name = annotate_dir+'/annotation_detail.csv'
	annotation_detail = open(annotation_detail_name, mode='w', newline='')
	annotation_writer = csv.writer(annotation_detail)
	annotation_writer.writerow(["seq_name", "seq_value", "seq_length",\
								"gene", "prokka_gene_name", "product", \
								"length", "start_pos", "end_pos", 'RGI_prediction_type'])
	ref_file_name = annotate_dir+ 'ref_neighborhood_sequences_1000_'+\
		datetime.datetime.now().strftime('%Y-%m-%d_%H-%M')+'.txt'
	ref_file = open(ref_file_name, 'w')
	ref_genome_files = extract_files('Experiments/1_1_1_limit/metagenome_data/', 'please provide the address of genome files')
	for genome_file in ref_genome_files:
		seq_list, _ = extract_amr_neighborhood_in_ref_genome(amr_seq, genome_file, 1000, 80)
		#AMR sequence might be in multiple places in ref genome so we have a list of
		#extracted neighborhood sequences instead of one sequence
		for index, seq in enumerate(seq_list):
			genome_name = os.path.splitext(os.path.basename(genome_file))[0]
			ref_file.write("> " + genome_name+'_'+str(index+1) + " :\n")
			ref_file.write(seq+"\n")
			seq_description = 'ref_'+genome_name+'_'+str(index+1)
			seq_info = annotate_sequence(seq+"\n", seq_description, annotate_dir+'/',
										'docker run -v `pwd`:/data staphb/prokka:latest ', True, False)
			for gene_info in seq_info:
					annotation_writer.writerow([seq_description, seq, len(seq),
							gene_info['gene'], gene_info['prokka_gene_name'],
							gene_info['product'], gene_info['length'],
							gene_info['start_pos'], gene_info['end_pos'],
							gene_info['RGI_prediction_type']])

	ref_file.close()
	annotation_detail.close()

def annotate_ref_genomes(amr_file, output_name, seq_length, ref_genome_files, amr_threshold,
						annotate_dir, prokka_prefix, use_RGI, RGI_include_loose,
						annotation_writer, trimmed_annotation_writer, gene_file, product_file):
	"""
	"""
	ref_up_info_list = []
	ref_down_info_list = []
	ref_amr_info_list = []
	#not necessarily unique values are stored here
	ref_seq_info_list = []
	amr_seq = retrieve_AMR(amr_file)
	ref_file_name = annotate_dir+'/ref_neighborhood_sequences'+output_name+'_'+\
		str(seq_length)+'_'+datetime.datetime.now().strftime('%Y-%m-%d_%H-%M')+'.txt'
	ref_file = open(ref_file_name, 'w')
	for genome_file in ref_genome_files:
		seq_list, _ = extract_amr_neighborhood_in_ref_genome(amr_seq, genome_file, seq_length, amr_threshold)
		#AMR sequence might be in multiple places in ref genome so we have a list of
		#extracted neighborhood sequences instead of one sequence
		for index, seq in enumerate(seq_list):
			genome_name = os.path.splitext(os.path.basename(genome_file))[0]
			ref_file.write("> " + genome_name+'_'+str(index+1) + " :\n")
			ref_file.write(seq+"\n")
			seq_description = 'ref_'+genome_name+'_'+str(index+1)
			seq_info = annotate_sequence(seq+"\n", seq_description, annotate_dir+'/',
										prokka_prefix, use_RGI, RGI_include_loose)
			myLine1 = myLine2 = seq_description + ':\t'
			found, ref_amr_info , ref_up_info, ref_down_info, seq_info = split_up_down_info(seq, seq_info)
			if found:
				if ref_up_info and not similar_seq_annotation_already_exist(ref_up_info, ref_up_info_list):
					ref_up_info_list.append(ref_up_info)
				if ref_down_info and not similar_seq_annotation_already_exist(ref_down_info, ref_down_info_list):
					ref_down_info_list.append(ref_down_info)
				ref_amr_info_list.append(ref_amr_info)
			else:
				logging.error("ERROR: no target amr was found in the ref sequence")
				import pdb; pdb.set_trace()
			for gene_info in seq_info:
				gene_info['seq_name'] = seq_description
				write_info_in_annotation_file(annotation_writer, trimmed_annotation_writer,
											gene_info, use_RGI, False)
				if gene_info['gene']=='':
					myLine1+='UNKNOWN---'
				else:
					myLine1+=gene_info['gene']+'---'
				myLine2+=gene_info['product']+'---'
			ref_seq_info_list.append(seq_info)
			gene_file.write(myLine1[:-3]+'\n')
			product_file.write(myLine2[:-3]+'\n')
	ref_file.close()
	logging.info("NOTE: The list of neighborhood sequences in reference genome(s) has\
	 		been stroed in " + ref_file_name)
	return ref_up_info_list, ref_down_info_list, ref_amr_info_list, ref_seq_info_list

def retrieve_neighborhood_seq(contig_name, up_limit, down_limit, ref_CAMI_file):
	"""
	"""
	#contig_name = contig_name.replace('_', '|')
	for record in SeqIO.parse(open(ref_CAMI_file,'r'),'fasta'):
		if contig_name in record.id:
			return str(record.seq[up_limit:min(down_limit, len(record.seq))])
	logging.error("ERROR: not able to find the contig "+contig_name)
	import pdb; pdb.set_trace()

def find_up_down_for_reverse_ref(gene_info, ref_seq_info, seq, ref_up_info_list,
						ref_down_info_list, index, up_limit, down_limit, use_RGI,
						gene_file, product_file, annotation_writer,
						trimmed_annotation_writer, to_be_removed, contig_index = 1):
	"""
	"""
	#upstream: find all annotations from the same contig that their end_pos <= down_limit
	ref_up_info = []
	counter = index+1
	for item in ref_seq_info[index+1:]:
		up_info = copy.deepcopy(item)
		if up_info['contig']!=gene_info['contig'] or max(up_info['start_pos'], up_info['end_pos']) > down_limit:
			break
		else:
			up_info['start_pos'] = len(seq)-(up_info['start_pos']-to_be_removed)+1
			up_info['end_pos'] = len(seq)-(up_info['end_pos']-to_be_removed)+1
			up_info['seq_value'] = seq
			ref_up_info.insert(0, up_info)
			counter+=1
	if ref_up_info and not similar_seq_annotation_already_exist(ref_up_info, ref_up_info_list):
		ref_up_info_list.append(ref_up_info)
	myLine1 = myLine2 = ''
	for gene_up_info in ref_up_info:
		seq_description = 'reference_'+gene_up_info['contig']+'_'+str(contig_index)
		gene_up_info['seq_name'] = seq_description
		write_info_in_annotation_file(annotation_writer, trimmed_annotation_writer,
								gene_up_info, use_RGI, False)
		if gene_up_info['gene']=='':
			myLine1+='UNKNOWN---'
		else:
			myLine1+=gene_up_info['gene']+'---'
		myLine2+=gene_up_info['product']+'---'
	gene_file.write(myLine1)
	product_file.write(myLine2)
	#amr gene
	seq_description = 'reference_'+gene_info['contig']+'_'+str(contig_index)
	gene_info['seq_name'] = seq_description
	gene_info['target_amr'] = 'yes'
	write_info_in_annotation_file(annotation_writer, trimmed_annotation_writer,
							gene_info, use_RGI, False)
	gene_file.write(gene_info['gene']+'---')
	product_file.write(gene_info['product']+'---')
	#downstream: find all annotations from the same contig that their start_pos >= up_limit
	ref_down_info = []
	for item in reversed(ref_seq_info[:index]):
		down_info = copy.deepcopy(item)
		if down_info['contig']!=gene_info['contig'] or min(down_info['start_pos'], down_info['end_pos']) < up_limit:
			break
		else:
			down_info['start_pos'] = len(seq)-(down_info['start_pos']-to_be_removed)+1
			down_info['end_pos'] = len(seq)-(down_info['end_pos']-to_be_removed)+1
			down_info['seq_value'] = seq
			ref_down_info.append(down_info)
	if ref_down_info and not similar_seq_annotation_already_exist(ref_down_info, ref_down_info_list):
		ref_down_info_list.append(ref_down_info)
	myLine1 = myLine2 = ''
	for gene_down_info in ref_down_info:
		seq_description = 'reference_'+gene_down_info['contig']+'_'+str(contig_index)
		gene_down_info['seq_name'] = seq_description
		write_info_in_annotation_file(annotation_writer, trimmed_annotation_writer,
								gene_down_info, use_RGI, False)
		if gene_down_info['gene']=='':
			myLine1+='UNKNOWN---'
		else:
			myLine1+=gene_down_info['gene']+'---'
		myLine2+=gene_down_info['product']+'---'
	gene_file.write(myLine1[:-3]+'\n')
	product_file.write(myLine2[:-3]+'\n')

	return ref_up_info_list, ref_down_info_list, ref_up_info+[gene_info]+ref_down_info, counter

def find_up_down_for_regular_ref(gene_info, ref_seq_info, seq, ref_up_info_list,
						ref_down_info_list, index, up_limit, down_limit, use_RGI,
						gene_file, product_file, annotation_writer,
						trimmed_annotation_writer, to_be_removed, contig_index = 1):
	"""
	"""
	#upstream: find all annotations from the same contig that their start_pos >= up_limit
	ref_up_info = []
	for item in reversed(ref_seq_info[:index]):
		up_info = copy.deepcopy(item)
		if up_info['contig']!=gene_info['contig'] or min(up_info['start_pos'], up_info['end_pos']) < up_limit:
			break
		else:
			up_info['start_pos'] -= to_be_removed
			up_info['end_pos'] -= to_be_removed
			up_info['seq_value'] = seq
			ref_up_info.insert(0, up_info)
	if ref_up_info and not similar_seq_annotation_already_exist(ref_up_info, ref_up_info_list):
		ref_up_info_list.append(ref_up_info)
	myLine1 = myLine2 = ''
	for gene_up_info in ref_up_info:
		seq_description = 'reference_'+gene_up_info['contig']+'_'+str(contig_index)
		gene_up_info['seq_name'] = seq_description
		write_info_in_annotation_file(annotation_writer, trimmed_annotation_writer,
								gene_up_info, use_RGI, False)
		if gene_up_info['gene']=='':
			myLine1+='UNKNOWN---'
		else:
			myLine1+=gene_up_info['gene']+'---'
		myLine2+=gene_up_info['product']+'---'
	gene_file.write(myLine1)
	product_file.write(myLine2)
	#amr gene
	seq_description = 'reference_'+gene_info['contig']+'_'+str(contig_index)
	gene_info['seq_name'] = seq_description
	gene_info['target_amr'] = 'yes'
	write_info_in_annotation_file(annotation_writer, trimmed_annotation_writer,
							gene_info, use_RGI, False)
	gene_file.write(gene_info['gene']+'---')
	product_file.write(gene_info['product']+'---')
	#downstream: find all annotations from the same contig that their end_pos <= down_limit
	ref_down_info = []
	counter = index+1
	for item in ref_seq_info[index+1:]:
		down_info = copy.deepcopy(item)
		if down_info['contig']!=gene_info['contig'] or max(down_info['start_pos'], down_info['end_pos']) > down_limit:
			break
		else:
			down_info['start_pos']-=to_be_removed
			down_info['end_pos']-=to_be_removed
			down_info['seq_value'] = seq
			ref_down_info.append(down_info)
			counter+=1
	if ref_down_info and not similar_seq_annotation_already_exist(ref_down_info, ref_down_info_list):
		ref_down_info_list.append(ref_down_info)
	myLine1 = myLine2 = ''
	for gene_down_info in ref_down_info:
		seq_description = 'reference_'+gene_down_info['contig']+'_'+str(contig_index)
		gene_down_info['seq_name'] = seq_description
		write_info_in_annotation_file(annotation_writer, trimmed_annotation_writer,
								gene_down_info, use_RGI, False)
		if gene_down_info['gene']=='':
			myLine1+='UNKNOWN---'
		else:
			myLine1+=gene_down_info['gene']+'---'
		myLine2+=gene_down_info['product']+'---'
	gene_file.write(myLine1[:-3]+'\n')
	product_file.write(myLine2[:-3]+'\n')

	return 	ref_up_info_list, ref_down_info_list, ref_up_info+[gene_info]+ref_down_info, counter

def extract_neighborhood_from_ref_annotation(amr_name, cami_info, seq_length, annotation_writer,
						trimmed_annotation_writer, gene_file, product_file, use_RGI):
	"""
	This is used for extraction annotation of neighborhoof of a given AMR
	from annotated gold standard fasta file available in CAMI.
	We look for the exact name as our target AMR; however, if nothing was found,
	we look for any gene from the same family of our target AMR!
	"""
	ref_up_info_list = []
	ref_down_info_list = []
	ref_amr_info_list = []
	ref_seq_info_list = []
	amr_family, ref_seq_info, ref_CAMI_file = cami_info
	index = 0
	contig_list =[]
	contig_index_list = []
	while index<len(ref_seq_info):
		gene_info = copy.deepcopy(ref_seq_info[index])
		reversed_direction = False
		gene = gene_info['gene'].strip().replace(' ','_')
		gene_name = ''.join(e for e in gene if e.isalpha() or e.isnumeric() or e=='_' or e=='-')
		if gene_name==amr_name and gene_info['RGI_prediction_type'] is not None:
			if gene_info['contig'] in contig_list:
				contig_index_list[contig_list.index(gene_info['contig'])]+=1
				contig_index = contig_index_list[contig_list.index(gene_info['contig'])]
			else:
				contig_list.append(gene_info['contig'])
				contig_index = 1
				contig_index_list.append(contig_index)
			if gene_info['start_pos'] > gene_info['end_pos']:
				reversed_direction = True
			start, end = min(gene_info['start_pos'], gene_info['end_pos']), max(gene_info['start_pos'], gene_info['end_pos'])
			up_limit = start-seq_length
			down_limit = end+seq_length
			seq = retrieve_neighborhood_seq(gene_info['contig'], max(up_limit, 0), down_limit, ref_CAMI_file)
			#!!!!!to_be_removed!!!!!!------seq_length--------S*******AMR******E
			to_be_removed = start - min(seq_length+1, start)
			#import pdb; pdb.set_trace()
			if (reversed_direction):
				seq=rc(seq)
				gene_info['start_pos'] = len(seq)-(gene_info['start_pos']-to_be_removed)+1
				gene_info['end_pos'] = len(seq)-(gene_info['end_pos']-to_be_removed)+1
				#gene_info['start_pos'], gene_info['end_pos'] = len(seq)-gene_info['end_pos']+1, len(seq)-gene_info['start_pos']+1
				gene_info['seq_value']=seq
				ref_amr_info_list.append(gene_info)
				ref_up_info_list, ref_down_info_list, new_ref_info, index = find_up_down_for_reverse_ref(
						gene_info, ref_seq_info, seq, ref_up_info_list, ref_down_info_list,
						index, up_limit, down_limit, use_RGI, gene_file, product_file,
						annotation_writer, trimmed_annotation_writer, to_be_removed, contig_index)
			else:
				gene_info['start_pos'] -= to_be_removed
				gene_info['end_pos'] -= to_be_removed
				gene_info['seq_value']=seq
				ref_amr_info_list.append(gene_info)
				ref_up_info_list, ref_down_info_list, new_ref_info, index = find_up_down_for_regular_ref(
						gene_info, ref_seq_info, seq, ref_up_info_list, ref_down_info_list,
						index, up_limit, down_limit, use_RGI, gene_file, product_file,
						annotation_writer, trimmed_annotation_writer, to_be_removed, contig_index)
			ref_seq_info_list.append(new_ref_info)
		else:
			index+=1
	if len(ref_amr_info_list)==0:
		ref_up_info_list, ref_down_info_list, ref_amr_info_list, ref_seq_info_list =\
			extract_neighborhood_from_ref_annotation_amr_family(amr_name, cami_info,
				seq_length, annotation_writer,trimmed_annotation_writer, gene_file, product_file, use_RGI)

	return ref_up_info_list, ref_down_info_list, ref_amr_info_list, ref_seq_info_list

def AMR_exists_in_rest_of_contig(amr_name, ref_seq_info, index):
	"""
	"""
	contig_name = ref_seq_info[index]['contig']
	for i, gene_info in enumerate(ref_seq_info[index+1:]):
		if gene_info['contig']==contig_name and gene_info['gene']==amr_name:
			return True, i+index+1
	return False, index

def extract_neighborhood_from_ref_annotation_amr_family(amr_name, cami_info, seq_length, annotation_writer,
						trimmed_annotation_writer, gene_file, product_file, use_RGI):
	"""
	This is used for extraction annotation of neighborhoof of a given AMR
	from annotated gold standard fasta file available in CAMI.
	Here, anywhere we can find a gene from the same family of our target AMR, we treat it as our AMR!
	"""
	ref_up_info_list = []
	ref_down_info_list = []
	ref_amr_info_list = []
	ref_seq_info_list = []
	amr_family, ref_seq_info, ref_CAMI_file = cami_info
	amr_family_list = [x.strip() for x in amr_family.split(';')]
	index = 0
	contig_list =[]
	contig_index_list = []
	while index<len(ref_seq_info):
		gene_info = copy.deepcopy(ref_seq_info[index])
		reversed_direction = False
		if gene_info['family'] is not None:
			gene_family_list = [x.strip() for x in gene_info['family'].split(';')]
			same_families = True
			for family in amr_family_list:
				if family not in gene_family_list:
					same_families = False
		else:
			same_families = False
		#if gene_info['family']==amr_family:
		if same_families:
			#first check if the gene name is different than amr_name and if that's the case
			#check the rest of this contig to see if you can find amr_name itself in which case
			# we directly move to that position in the contig and ignore other potential
			#cases from the same family but with a different name
			#This is to avoid this problem:
			# A --- B (AMR) --- C
			#ref_CAMI_annotated:   A (AMR) --- B --- C
			#because A nad B are from the same family, while we are looking for B as the target AMR
			if gene_info['gene']!=amr_name:
				found, index = AMR_exists_in_rest_of_contig(amr_name, ref_seq_info, index)
				if found:
					gene_info = copy.deepcopy(ref_seq_info[index])
			if gene_info['contig'] in contig_list:
				contig_index_list[contig_list.index(gene_info['contig'])]+=1
				contig_index = contig_index_list[contig_list.index(gene_info['contig'])]
			else:
				contig_list.append(gene_info['contig'])
				contig_index = 1
				contig_index_list.append(contig_index)
			if gene_info['start_pos'] > gene_info['end_pos']:
				reversed_direction = True
			start, end = min(gene_info['start_pos'], gene_info['end_pos']), max(gene_info['start_pos'], gene_info['end_pos'])
			up_limit = start-seq_length
			down_limit = end+seq_length
			seq = retrieve_neighborhood_seq(gene_info['contig'], max(up_limit, 0), down_limit, ref_CAMI_file)
			#!!!!!to_be_removed!!!!!!------seq_length--------S*******AMR******E
			to_be_removed = start - min(seq_length+1, start)
			#import pdb; pdb.set_trace()
			if (reversed_direction):
				seq=rc(seq)
				gene_info['start_pos'] = len(seq)-(gene_info['start_pos']-to_be_removed)+1
				gene_info['end_pos'] = len(seq)-(gene_info['end_pos']-to_be_removed)+1
				#gene_info['start_pos'], gene_info['end_pos'] = len(seq)-gene_info['end_pos']+1, len(seq)-gene_info['start_pos']+1
				gene_info['seq_value']=seq
				ref_amr_info_list.append(gene_info)
				ref_up_info_list, ref_down_info_list, new_ref_info, index = find_up_down_for_reverse_ref(
						gene_info, ref_seq_info, seq, ref_up_info_list, ref_down_info_list,
						index, up_limit, down_limit, use_RGI, gene_file, product_file,
						annotation_writer, trimmed_annotation_writer, to_be_removed, contig_index)
			else:
				gene_info['start_pos'] -= to_be_removed
				gene_info['end_pos'] -= to_be_removed
				gene_info['seq_value']=seq
				ref_amr_info_list.append(gene_info)
				ref_up_info_list, ref_down_info_list, new_ref_info, index = find_up_down_for_regular_ref(
						gene_info, ref_seq_info, seq, ref_up_info_list, ref_down_info_list,
						index, up_limit, down_limit, use_RGI, gene_file, product_file,
						annotation_writer, trimmed_annotation_writer, to_be_removed, contig_index)
			ref_seq_info_list.append(new_ref_info)
		else:
			index+=1

	return ref_up_info_list, ref_down_info_list, ref_amr_info_list, ref_seq_info_list

def extract_ref_neighborhood_and_annotate(amr_file, amr_name, cami_info, seq_length,
								amr_threshold, annotate_dir, prokka_prefix, use_RGI,
								RGI_include_loose, annotation_writer, trimmed_annotation_writer,
								gene_file, product_file):
	"""
	AMR alignment --> neighborhood extraction --> annotation for CAMI gold standard contig list
	"""
	_, _, ref_CAMI_file = cami_info
	ref_up_info_list = []
	ref_down_info_list = []
	ref_amr_info_list = []
	ref_seq_info_list = []
	amr_seq = retrieve_AMR(amr_file)
	ref_file_name = annotate_dir+'/ref_neighborhood_sequences_'+amr_name+'_'+\
		str(seq_length)+'_'+datetime.datetime.now().strftime('%Y-%m-%d_%H-%M')+'.txt'
	ref_file = open(ref_file_name, 'w')
	seq_list, contig_name_list = extract_amr_neighborhood_in_ref_genome(amr_seq, ref_CAMI_file, seq_length, amr_threshold)
	#AMR sequence might be in multiple places in ref genome so we have a list of
	#extracted neighborhood sequences instead of one sequence
	unique_contig_name_list = []
	contig_index_list = []
	for index, seq in enumerate(seq_list):
		contig_name = contig_name_list[index]
		if contig_name in unique_contig_name_list:
			myindex = unique_contig_name_list.index(contig_name)
			contig_index_list[myindex]+=1
			contig_name+='_'+str(contig_index_list[myindex])
		else:
			unique_contig_name_list.append(contig_name)
			contig_index_list.append(1)
			contig_name+='_1'
		ref_file.write("> " + contig_name+" :\n")
		ref_file.write(seq+"\n")
		seq_description = 'ref_'+contig_name
		seq_info = annotate_sequence(seq+"\n", seq_description, annotate_dir+'/',
									prokka_prefix, use_RGI, RGI_include_loose)
		myLine1 = myLine2 = seq_description + ':\t'
		found, ref_amr_info , ref_up_info, ref_down_info, seq_info = split_up_down_info(seq, seq_info)
		if found:
			if ref_up_info and not similar_seq_annotation_already_exist(ref_up_info, ref_up_info_list):
				ref_up_info_list.append(ref_up_info)
			if ref_down_info and not similar_seq_annotation_already_exist(ref_down_info, ref_down_info_list):
				ref_down_info_list.append(ref_down_info)
			ref_amr_info_list.append(ref_amr_info)
		else:
			logging.error("ERROR: no target amr was found in the ref sequence")
			return [], [], []
			#import pdb; pdb.set_trace()
		for gene_info in seq_info:
			gene_info['seq_name'] = seq_description
			write_info_in_annotation_file(annotation_writer, trimmed_annotation_writer,
										gene_info, use_RGI, False)
			if gene_info['gene']=='':
				myLine1+='UNKNOWN---'
			else:
				myLine1+=gene_info['gene']+'---'
			myLine2+=gene_info['product']+'---'
		ref_seq_info_list.append(seq_info)
		gene_file.write(myLine1[:-3]+'\n')
		product_file.write(myLine2[:-3]+'\n')
	ref_file.close()
	logging.info("NOTE: The list of neighborhood sequences in reference genome(s) has\
	 		been stroed in " + ref_file_name)
	return ref_up_info_list, ref_down_info_list, ref_amr_info_list, ref_seq_info_list

def extract_graph_seqs_annotation_parallel(amr_file, path_info_file, neighborhood_seq_file,
					annotate_dir, core_num, ref_genomes_available, prokka_prefix,
					use_RGI, RGI_include_loose,ref_up_info_list, ref_down_info_list,
					annotation_writer, trimmed_annotation_writer, gene_file, product_file):
	"""
	"""
	#Read path_info from the file
	path_info_list = []
	if path_info_file!=-1:
		path_info_list = read_path_info_file(path_info_file)
	#find the list of all extracted sequences
	logging.info('Reading '+ neighborhood_seq_file + ' for '+ amr_file)
	sequence_list = []
	counter = 1
	with open(neighborhood_seq_file, 'r') as read_obj:
		for line in read_obj:
			if line.startswith('>') or line.startswith('Path') or line.startswith('The'):
				continue
			sequence_list.append((counter, line))
			counter+=1
	#Parallel annotation
	p_annotation = partial(extract_seq_annotation, annotate_dir, prokka_prefix,
							use_RGI, RGI_include_loose)
	with Pool(core_num) as p:
		seq_info_list = p.map(p_annotation, sequence_list)
	#Further processing of result of parallel annotation
	all_seq_info_list =[]
	for i, seq_pair in enumerate(sequence_list):
		counter, line = seq_pair
		seq_description = 'extracted'+str(counter)
		seq_info = seq_info_list[i]
		#extract amr from seq_info
		amr_found, amr_info, up_info, down_info, seq_info = split_up_down_info(line[:-1], seq_info)
		if not amr_found:
			logging.error("ERROR: no target amr was found in the extracted sequence")
			import pdb; pdb.set_trace()
		#calculate coverage for the genes available in the annotation
		if path_info_list:
			coverage_list = find_gene_coverage(seq_info, path_info_list[counter-1])
		#Check if this annotation has already been found
		found = seq_annotation_already_exist(seq_info, all_seq_info_list)
		#If it's a novel sequence correct annotation (if applicable) for cases that RGI doesn't have a hit but Prokka has
		if not found:
			all_seq_info_list.append(seq_info)
			if ref_genomes_available:
				seq_info = add_RGI_loose(line[:-1], seq_info, ref_up_info_list, ref_down_info_list)
		myLine1 = myLine2 = seq_description +':\t'
		#write annotation onfo into the files
		for j, gene_info in enumerate(seq_info):
			coverage = coverage_list[j]
			gene_info['coverage'] = coverage
			gene_info['seq_name'] = seq_description
			write_info_in_annotation_file(annotation_writer, trimmed_annotation_writer,
										gene_info, use_RGI, found)
			if gene_info['gene']=='':
				myLine1+='UNKNOWN---'
			else:
				myLine1+=gene_info['gene']+'---'
			myLine2+=gene_info['product']+'---'
		gene_file.write(myLine1[:-3]+'\n')
		product_file.write(myLine2[:-3]+'\n')
	return all_seq_info_list

def neighborhood_annotation_parallel(amr_file, ref_genome_files, neighborhood_seq_file, \
								path_info_file, seq_length,\
								output_dir, prokka_prefix, use_RGI = True,\
								RGI_include_loose = False, output_name ='',
								amr_threshold = 95, ref_genomes_available = True,
								core_num = 4, cami_info = None):
	"""
	To annotate reference genomes (a piece extracted around the AMR gene) as well as
		extracted neighborhood sequences from assembly graph, summarize the results
		in a couple of formats and visualize them
	Parameters:
		amr_file:	the address of the file containing the AMR sequence
		ref_genome_files:	the list of address of reference genomes after
			inserting the AMR gene in them
		neighborhood_seq_file:	the address of the file containing all extracted
		 	neighborhood sequences from assembly graph
		seq_length:	the length of neighborhood sequence extracted from each side of
			the AMR sequence (downstream and up stream)
		output_dir:	the path for the output directory
		prokka_prefix: it's used to run prokka properly via docker or conda
		use_RGI:	RGI annotations incorporated for AMR annotation
		RGI_include_loose: Whether to include loose annotaions in RGI
		output_name:the name used to distinguish different output files usually based on the name of AMR
		amr_threshold: the threshold used for identity and coverage
	Return:
		the address of files stroing annotation information (annotation_detail_name,
			trimmed_annotation_info, gene_file_name, product_file_name, visual_annotation)
	"""
	logging.debug('Started annotation for '+output_name)
	# initializing required files and directories
	annotate_dir = output_dir+ANNOTATION_DIR+'/'+ANNOTATION_DIR+'_'+str(seq_length)+'/annotation'+output_name+'_'+str(seq_length)
	if os.path.exists(annotate_dir):
		try:
			shutil.rmtree(annotate_dir)
		except OSError as e:
			logging.error("Error: %s - %s." % (e.filename, e.strerror))
			#print("Error: %s - %s." % (e.filename, e.strerror))
	os.makedirs(annotate_dir)
	annotation_detail_name = annotate_dir+'/annotation_detail'+output_name+'.csv'
	trimmed_annotation_info_name = annotate_dir+'/trimmed_annotation_info'+output_name+'.csv'
	annotation_detail = open(annotation_detail_name, mode='w', newline='')
	trimmed_annotation_info = open(trimmed_annotation_info_name, mode='w', newline='')
	annotation_writer = csv.writer(annotation_detail)
	trimmed_annotation_writer = csv.writer(trimmed_annotation_info)
	gene_info = {'seq_value':'seq_value', 'gene':'gene', 'prokka_gene_name':'prokka_gene_name',
				'product':'product', 'length':'length', 'start_pos':'start_pos',
				'end_pos':'end_pos', 'RGI_prediction_type':'RGI_prediction_type',
				'coverage':'coverage', 'family':'family', 'seq_name':'seq_name',
				'target_amr':'target_amr'}
	write_info_in_annotation_file(annotation_writer, trimmed_annotation_writer,
								gene_info, use_RGI, False, 'seq_length')
	gene_file_name = annotate_dir+'/seq_comparison_genes'+output_name+'.txt'
	gene_file = open(gene_file_name, 'w')
	product_file_name = annotate_dir+'/seq_comparison_products'+output_name+'.txt'
	product_file = open(product_file_name, 'w')

	#find neighborhood sequence in ref genome(s) and annotate using prokka and RGI
	ref_up_info_list = []
	ref_down_info_list = []
	ref_amr_info_list = []
	ref_seq_info_list = []
	if (ref_genomes_available):
		ref_up_info_list, ref_down_info_list, ref_amr_info_list, ref_seq_info_list = annotate_ref_genomes(
					amr_file, output_name, seq_length, ref_genome_files, amr_threshold,
					annotate_dir, prokka_prefix, use_RGI, RGI_include_loose,
					annotation_writer, trimmed_annotation_writer, gene_file, product_file)
	elif cami_info:
		if CAMI_REF_NEIGHBORHOOD_METHOD == 1:
			ref_up_info_list, ref_down_info_list, ref_amr_info_list, ref_seq_info_list =\
			 					extract_neighborhood_from_ref_annotation_amr_family(
								output_name[1:], cami_info, seq_length, annotation_writer,
								trimmed_annotation_writer, gene_file, product_file, use_RGI)
		elif CAMI_REF_NEIGHBORHOOD_METHOD == 2:
			ref_up_info_list, ref_down_info_list, ref_amr_info_list, ref_seq_info_list =\
			 					extract_neighborhood_from_ref_annotation(
								output_name[1:], cami_info, seq_length, annotation_writer,
								trimmed_annotation_writer, gene_file, product_file, use_RGI)
		elif CAMI_REF_NEIGHBORHOOD_METHOD == 3:
			ref_up_info_list, ref_down_info_list, ref_amr_info_list, ref_seq_info_list =\
								extract_ref_neighborhood_and_annotate(
								amr_file, output_name[1:], cami_info, seq_length,
								amr_threshold, annotate_dir, prokka_prefix, use_RGI,
								RGI_include_loose, annotation_writer, trimmed_annotation_writer,
								gene_file, product_file)

	#annotate the sequences extraced from assembly graph
	all_seq_info_list = extract_graph_seqs_annotation_parallel(amr_file, path_info_file, neighborhood_seq_file,
									annotate_dir, core_num, ref_genomes_available,
									prokka_prefix, use_RGI, RGI_include_loose,
									ref_up_info_list, ref_down_info_list,
									annotation_writer, trimmed_annotation_writer,
									gene_file, product_file)
	logging.info("NOTE: The comparison of neighborhood sequences are available in " +\
	 		annotation_detail_name+", "+gene_file_name+", "+product_file_name)
	annotation_detail.close()
	trimmed_annotation_info.close()
	gene_file.close()
	product_file.close()

	return all_seq_info_list, ref_seq_info_list,ref_up_info_list, ref_down_info_list,\
			ref_amr_info_list, trimmed_annotation_info_name

def extract_graph_seqs_annotation(amr_file, path_info_file, neighborhood_seq_file,
					annotate_dir, ref_genomes_available, prokka_prefix,
					use_RGI, RGI_include_loose,ref_up_info_list, ref_down_info_list,
					annotation_writer, trimmed_annotation_writer, gene_file, product_file):
	"""
	"""
	counter = 1
	all_seq_info_list =[]
	#Read path_info from the file
	path_info_list = []
	if path_info_file!=-1:
		path_info_list = read_path_info_file(path_info_file)
	logging.info('Reading '+ neighborhood_seq_file + ' for '+ amr_file)
	#print('Reading '+ neighborhood_seq_file + ' for '+ amr_file )
	with open(neighborhood_seq_file, 'r') as read_obj:
		for line in read_obj:
			if line.startswith('>') or line.startswith('Path') or line.startswith('The'):
				continue
			seq_description = 'extracted'+str(counter)
			seq_info = annotate_sequence(line, seq_description, annotate_dir+'/',
											prokka_prefix, use_RGI, RGI_include_loose)
			#extract amr from seq_info
			amr_found, amr_info, up_info, down_info, seq_info = split_up_down_info(line[:-1], seq_info)
			if not amr_found:
				logging.error("ERROR: no target amr was found in the extracted sequence")
				import pdb; pdb.set_trace()
			#calculate the coverage of annotated genes
			if path_info_list:
				coverage_list = find_gene_coverage(seq_info, path_info_list[counter-1])
			#Check if this annotation has already been found
			found = seq_annotation_already_exist(seq_info, all_seq_info_list)
			#If it's a novel sequence correct annotation if possible
			if not found:
				all_seq_info_list.append(seq_info)
				if ref_genomes_available:
					seq_info = add_RGI_loose(line[:-1], seq_info, ref_up_info_list, ref_down_info_list)
			myLine1 = myLine2 = seq_description +':\t'
			#write annotation onfo into the files
			for j, gene_info in enumerate(seq_info):
				coverage = coverage_list[j]
				gene_info['coverage'] = coverage
				gene_info['seq_name'] = seq_description
				write_info_in_annotation_file(annotation_writer, trimmed_annotation_writer,
											gene_info, use_RGI, found)
				if gene_info['gene']=='':
					myLine1+='UNKNOWN---'
				else:
					myLine1+=gene_info['gene']+'---'
				myLine2+=gene_info['product']+'---'
			gene_file.write(myLine1[:-3]+'\n')
			product_file.write(myLine2[:-3]+'\n')
			counter+=1
	return all_seq_info_list

def neighborhood_annotation(amr_file, ref_genome_files, neighborhood_seq_file, \
								path_info_file, seq_length,\
								output_dir, prokka_prefix, use_RGI = True,\
								RGI_include_loose = False, output_name ='',
								amr_threshold = 95, ref_genomes_available = True,
								cami_info = None):
	"""
	To annotate reference genomes (a piece extracted around the AMR gene) as well as
		extracted neighborhood sequences from assembly graph, summarize the results
		in a couple of formats and visualize them.
	Parameters:
		amr_file:	the address of the file containing the AMR sequence
		ref_genome_files:	the list of address of reference genomes after
			inserting the AMR gene in them
		neighborhood_seq_file:	the address of the file containing all extracted
		 	neighborhood sequences from assembly graph
		seq_length:	the length of neighborhood sequence extracted from each side of
			the AMR sequence (downstream and up stream)
		output_dir:	the path for the output directory
		prokka_prefix: it's used to run prokka properly via docker or conda
		use_RGI:	RGI annotations incorporated for AMR annotation
		RGI_include_loose: Whether to include loose annotaions in RGI
		output_name:the name used to distinguish different output files usually based on the name of AMR
		amr_threshold: the threshold used for identity and coverage
	Return:
		the address of files storing annotation information (annotation_detail_name,
			trimmed_annotation_info, gene_file_name, product_file_name, visual_annotation)
	"""
	logging.debug('Started annotation for '+output_name)
	# initializing required files and directories
	annotate_dir = output_dir+ANNOTATION_DIR+'/'+ANNOTATION_DIR+'_'+str(seq_length)+'/annotation'+output_name+'_'+str(seq_length)
	if os.path.exists(annotate_dir):
		try:
			shutil.rmtree(annotate_dir)
		except OSError as e:
			logging.error("Error: %s - %s." % (e.filename, e.strerror))
	os.makedirs(annotate_dir)
	annotation_detail_name = annotate_dir+'/annotation_detail'+output_name+'.csv'
	trimmed_annotation_info_name = annotate_dir+'/trimmed_annotation_info'+output_name+'.csv'
	annotation_detail = open(annotation_detail_name, mode='w', newline='')
	trimmed_annotation_info = open(trimmed_annotation_info_name, mode='w', newline='')
	#annotation_writer = csv.writer(annotation_detail, delimiter=',', quoting=csv.QUOTE_MINIMAL)
	annotation_writer = csv.writer(annotation_detail)
	trimmed_annotation_writer = csv.writer(trimmed_annotation_info)
	gene_info = {'seq_value':'seq_value', 'gene':'gene', 'prokka_gene_name':'prokka_gene_name',
				'product':'product', 'length':'length', 'start_pos':'start_pos',
				'end_pos':'end_pos', 'RGI_prediction_type':'RGI_prediction_type',
				'coverage':'coverage', 'family':'family', 'seq_name':'seq_name',
				'target_amr':'target_amr'}
	write_info_in_annotation_file(annotation_writer, trimmed_annotation_writer,
								gene_info, use_RGI, False, 'seq_length')
	gene_file_name = annotate_dir+'/seq_comparison_genes'+output_name+'.txt'
	gene_file = open(gene_file_name, 'w')
	product_file_name = annotate_dir+'/seq_comparison_products'+output_name+'.txt'
	product_file = open(product_file_name, 'w')
	#find neighborhood sequence in ref genome(s) and annotate using prokka and RGI
	if (ref_genomes_available):
		ref_up_info_list, ref_down_info_list, ref_amr_info_list, ref_seq_info_list = annotate_ref_genomes(
					amr_file, output_name, seq_length, ref_genome_files, amr_threshold,
					annotate_dir, prokka_prefix, use_RGI, RGI_include_loose,
					annotation_writer, trimmed_annotation_writer, gene_file, product_file)
	elif cami_info:
		if CAMI_REF_NEIGHBORHOOD_METHOD == 1:
			ref_up_info_list, ref_down_info_list, ref_amr_info_list, ref_seq_info_list =\
								extract_neighborhood_from_ref_annotation_amr_family(
								output_name[1:], cami_info, seq_length, annotation_writer,
								trimmed_annotation_writer, gene_file, product_file, use_RGI)
		elif CAMI_REF_NEIGHBORHOOD_METHOD == 2:
			ref_up_info_list, ref_down_info_list, ref_amr_info_list, ref_seq_info_list =\
								extract_neighborhood_from_ref_annotation(
								output_name[1:], cami_info, seq_length, annotation_writer,
								trimmed_annotation_writer, gene_file, product_file, use_RGI)
		elif CAMI_REF_NEIGHBORHOOD_METHOD == 3:
			ref_up_info_list, ref_down_info_list, ref_amr_info_list, ref_seq_info_list =\
								extract_ref_neighborhood_and_annotate(
								amr_file, output_name[1:], cami_info, seq_length,
								amr_threshold, annotate_dir, prokka_prefix, use_RGI,
								RGI_include_loose, annotation_writer, trimmed_annotation_writer,
								gene_file, product_file)

	#annotate the sequences extraced from assembly graph
	all_seq_info_list = extract_graph_seqs_annotation(amr_file, path_info_file, neighborhood_seq_file,
					annotate_dir, ref_genomes_available, prokka_prefix,
					use_RGI, RGI_include_loose,ref_up_info_list, ref_down_info_list,
					annotation_writer, trimmed_annotation_writer, gene_file, product_file)
	annotation_detail.close()
	trimmed_annotation_info.close()
	gene_file.close()
	product_file.close()
	logging.info("NOTE: The comparison of neighborhood sequences are available in " +\
	 		annotation_detail_name+", "+gene_file_name+", "+product_file_name)

	return all_seq_info_list, ref_seq_info_list,ref_up_info_list, ref_down_info_list,\
			ref_amr_info_list, trimmed_annotation_info_name

def create_fasta_file(seq, output_dir, comment = "> sequence:\n", file_name = 'temp'):
	"""
	To create a fasta file for a sequence
	Parameters:
		seq: the sequence to be written into the file
		output_dir: the output directory address
		comment: the comment to be written into fasta file
		file_name: the name of the fasta file
	Return:
		the address of the fasta file
	"""
	myfile_name = output_dir+file_name+'.fasta'
	if os.path.isfile(myfile_name):
		os.remove(myfile_name)
	myfile = open(myfile_name, 'w')
	myfile.write(comment)
	myfile.write(seq)
	myfile.close()
	return  myfile_name


def is_there_amr_in_graph(amr_file, amr_name, gfa_file, output_dir, bandage_path, threshold =  99):
	"""
	To call bandage+blast and check if the amr sequence can be found in the assembly graph
	Parameters:
		amr_file: the address of the query file
		amr_name: the name of the AMR sequence
		gfa_file: the address of the assembly graph
		output_dir: the address of the output directory
		bandage_path: the address of bandage executable file
		threshold: the threshold for coverage and identity
	Return:
		a boolean value which is True if amr_file was found in gfa_file and the
		list of paths returned by bandage+blast in which the coverage and identiry
		are greater/equal than/to threshold
	"""
	output_name=output_dir+amr_name+'_align_'+datetime.datetime.now().strftime('%Y-%m-%d_%H-%M-%S')
	if os.path.isfile(output_name+'.tsv'):
		os.remove(output_name+'.tsv')
	command = bandage_path +' querypaths '+gfa_file+' '+amr_file+' '+output_name
	os.system(command)

	paths_info = []
	found = False
	with open(output_name+".tsv") as tsvfile:
		reader = csv.reader(tsvfile, delimiter='\t')
		#skip the header
		next(reader)
		for row in reader:
			coverage = float(re.sub('[%]','',row[3]))
			identity = float(re.sub('[%]','',row[5]))
			if int(coverage) >= threshold and int(identity)>=threshold:
				found = True
				cell_info = row[1].strip()
				nodes, orientation_list, start_pos, end_pos = extract_nodes_in_path(cell_info)
				path_info = {'nodes':nodes, 'orientations':orientation_list,
								'start_pos':start_pos, 'end_pos':end_pos}
				paths_info.append(path_info)
	if not found:
		os.remove(output_name+'.tsv')
		logging.debug(amr_name+' not found in the graph!')
	else:
		logging.debug(amr_name+' found!!!')
	return found, paths_info, output_name+".tsv"

def amr_path_overlap(found_amr_paths, new_paths, new_amr_len, overlap_percent = 95):
	"""
	To check if all paths found for the new AMR seq overlap significantly (greater/equal
	 than/to overlap percent) with the already found paths for other AMRs
	Parameters:
	 	found_amr_paths:  	the paths already found for AMR genes
		new_paths: 			the paths found for the new AMR gene
		overlap_percent:	the threshold for overlap
	Return:
		False only if every paths in new_paths have overlap with at least one path in found_amr_paths
		True if we can find at least one path that is unique and not available in found_amr_paths
		Also, it returns the list of indeces from found_amr_paths that had overlap with a path in new_paths
	"""
	id_list = []
	for new_path in new_paths:
		found = False
		for i, paths in enumerate(found_amr_paths):
			for path in paths:
				if path['nodes'] == new_path['nodes'] and path['orientations']==new_path['orientations']:
					# for now we just check overlaps when they are in the same node(s)
					diff_length = max(path['start_pos']-new_path['start_pos'],0)+max(new_path['end_pos']-path['end_pos'],0)
					percent = (1 - (float(diff_length)/(new_amr_len)))*100
					if percent >= overlap_percent:
						found = True
						id_list.append(i)
						break
			if found:
				break
	if len(id_list)==len(new_paths):
		return True, id_list
	return False, None

def extract_amr_object_and_find_parallel(gfa_file, output_dir, align_dir, bandage_path,
											amr_threshold, amr_object):
	"""
	"""
	amr_seq, amr_name = amr_object
	logging.info('Checking if AMR "'+amr_name+'" exists in the assembly graph...')
	#create a fasta file for it
	amr_file = create_fasta_file(amr_seq, output_dir, file_name=amr_name)
	#Run Bandage+BLAST
	found, paths_info, tsv_file = is_there_amr_in_graph(amr_file, amr_name, gfa_file,
									align_dir, bandage_path, amr_threshold)
	#Remove temporary AMR file
	if os.path.isfile(amr_file):
		os.remove(amr_file)
	return found, paths_info, tsv_file

def find_all_amr_in_graph_parallel(gfa_file, output_dir, amr_sequences_file,
									bandage_path, amr_threshold, core_num):
	"""
	To go over a list of AMR sequences (amr_sequences_file) and run bandage+blast
	to check if any of them exists in the assembly graph (gfa_file)
	Parameters:
		gfa_file: the address of the assembly graph
		output_dir: the address of the output directory
		amr_sequences_file: the address of the file containing the sequence of all AMRs from CARD
		bandage_path: the address of bandage executable file
		amr_threshold: the threshold for coverage and identity
	Return:

	"""
	align_dir = output_dir+AMR_DIR_NAME+AMR_ALIGN_DIR
	if not os.path.exists(align_dir):
		os.makedirs(align_dir)

	amr_name = ''
	found_amr_seqs = []
	found_amr_names = []
	found_amr_paths = []
	#Read AMR sequences one by one
	amr_objects_parallel = []
	amr_objects = []
	with open(amr_sequences_file) as fp:
		for line in fp:
			if line.startswith('>'):
				amr_name = line
				continue
			amr_name_processed1 = amr_name.split('[')[0].split('|')[-1].strip().replace(' ','_')
			amr_name_processed = ''.join(e for e in amr_name_processed1 if e.isalpha() or e.isnumeric() or e=='_' or e=='-')
			amr_object={'seq':line, 'name':amr_name, 'p_name1':amr_name_processed1}
			amr_objects_parallel.append((line, amr_name_processed))
			amr_objects.append(amr_object)
	#parallel run Bandage+BLAST
	p_find_amr = partial(extract_amr_object_and_find_parallel, gfa_file, output_dir,
						align_dir, bandage_path, amr_threshold)
	with Pool(core_num) as p:
		lists = p.map(p_find_amr, amr_objects_parallel)
	found, paths_info, tsv_files = zip(*lists)

	align_files = []
	#process the result of parallel processes
	for i in range(len(amr_objects)):
		if found[i]:
			logging.debug(amr_objects[i]['name']+'('+amr_objects[i]['p_name1']+' or '+ amr_objects_parallel[i][1]+') was found: '+tsv_files[i])
			overlap, amr_ids =  amr_path_overlap(found_amr_paths, paths_info[i],
									len(amr_objects[i]['seq'])-1)
			if not overlap:
				logging.debug('No overlap for '+amr_objects_parallel[i][1])
				found_amr_seqs.append(amr_objects[i]['seq'])
				amr_info = {'name':amr_objects[i]['name'], 'overlap_list':[]}
				found_amr_names.append(amr_info)
				found_amr_paths.append(paths_info[i])
				align_files.append(tsv_files[i])
			else:
				# add this AMR to the right group of AMRs all having overlaps
				for id in amr_ids:
					if amr_objects[i]['p_name1'] not in found_amr_names[id]['overlap_list']:
						found_amr_names[id]['overlap_list'].append(amr_objects[i]['p_name1'])

	# write information (the sequence of found AMRs that don't have overlaped paths with others
	# + the list of groups in each all AMRs have overlaped paths) into files
	AMR_dir = output_dir+AMR_DIR_NAME+AMR_SEQ_DIR
	if not os.path.exists(AMR_dir):
		os.makedirs(AMR_dir)
	overlap_file_name = output_dir+AMR_DIR_NAME+AMR_OVERLAP_FILE
	overlap_file = open(overlap_file_name, 'w')
	amr_files = []
	for i, seq in enumerate(found_amr_seqs):
		amr_name1 = found_amr_names[i]['name'].split('[')[0].split('|')[-1].strip().replace(' ','_')
		amr_name = ''.join(e for e in amr_name1 if e.isalpha() or e.isnumeric() or e=='_' or e=='-')
		amr_file = create_fasta_file(seq, AMR_dir, found_amr_names[i]['name'], amr_name)
		amr_files.append(amr_file)
		if found_amr_names[i]['overlap_list']:
			overlap_file.write(amr_name1+":\n")
			overlap_file.write(', '.join(e for e in found_amr_names[i]['overlap_list']))
			overlap_file.write("\n")
	overlap_file.close()

	return amr_files, align_files

def find_all_amr_in_graph(gfa_file, output_dir, amr_sequences_file, bandage_path, amr_threshold):
	"""
	To go over a list of AMR sequences (amr_sequences_file) and run bandage+blast
	to check if any of them exists in the assembly graph (gfa_file)
	Parameters:
		gfa_file: the address of the assembly graph
		output_dir: the address of the output directory
		amr_sequences_file: the address of the file containing the sequence of all AMRs from CARD
		bandage_path: the address of bandage executable file
		amr_threshold: the threshold for coverage and identity
	Return:

	"""
	align_dir = output_dir+AMR_DIR_NAME+AMR_ALIGN_DIR
	if not os.path.exists(align_dir):
		os.makedirs(align_dir)

	amr_name = ''
	found_amr_seqs = []
	found_amr_names = []
	found_amr_paths = []
	align_files = []
	#Read AMR sequences one by one
	with open(amr_sequences_file) as fp:
		for line in fp:
			if line.startswith('>'):
				amr_name = line
				continue
			#create a fasta file for it
			amr_file = create_fasta_file(line, output_dir)
			#Run Bandage+BLAST
			amr_name_processed1 = amr_name.split('[')[0].split('|')[-1].strip().replace(' ','_')
			amr_name_processed = ''.join(e for e in amr_name_processed1 if e.isalpha() or e.isnumeric() or e=='_' or e=='-')
			logging.info('Checking if AMR "'+amr_name_processed+'" exists in the assembly graph...')
			found, paths_info, tsv_file = is_there_amr_in_graph(amr_file, amr_name_processed, gfa_file,
											align_dir, bandage_path, amr_threshold)
			if found:
				overlap, amr_ids =  amr_path_overlap(found_amr_paths, paths_info,
										len(line)-1)
				if not overlap:
					found_amr_seqs.append(line)
					amr_info = {'name':amr_name, 'overlap_list':[]}
					found_amr_names.append(amr_info)
					found_amr_paths.append(paths_info)
					align_files.append(tsv_file)
				else:
					# add this AMR to the right group of AMRs all having overlaps
					for id in amr_ids:
						if amr_name_processed1 not in found_amr_names[id]['overlap_list']:
							found_amr_names[id]['overlap_list'].append(amr_name_processed1)
	import pdb; pdb.set_trace()
	# write information (the sequence of found AMRs that don't have overlaped paths with others
	# + the list of groups in each all AMRs have overlaped paths) into files
	AMR_dir = output_dir+AMR_DIR_NAME+AMR_SEQ_DIR
	if not os.path.exists(AMR_dir):
		os.makedirs(AMR_dir)
	overlap_file_name = output_dir+AMR_DIR_NAME+AMR_OVERLAP_FILE
	overlap_file = open(overlap_file_name, 'w')
	amr_files = []
	for i, seq in enumerate(found_amr_seqs):
		amr_name1 = found_amr_names[i]['name'].split('[')[0].split('|')[-1].strip().replace(' ','_')
		amr_name = ''.join(e for e in amr_name1 if e.isalpha() or e.isnumeric() or e=='_' or e=='-')
		amr_file = create_fasta_file(seq, AMR_dir, found_amr_names[i]['name'], amr_name)
		amr_files.append(amr_file)
		if found_amr_names[i]['overlap_list']:
			overlap_file.write(amr_name1+":\n")
			overlap_file.write(', '.join(e for e in found_amr_names[i]['overlap_list']))
			overlap_file.write("\n")
	overlap_file.close()

	return amr_files, align_files

def find_corrsponding_seq_path_file(amr_name, sequences_file_names, path_info_file_names, seq_length):
	"""
	To return the name of a sequence file (from sequences_file_names)
	and the name of a path file (from path_info_file_names)
	dedicated to sequences extracted with a given length (seq_length)
	from a given amr sequence (amr_name)
	"""
	seq_file = -1
	for file_name in sequences_file_names:
		if SEQ_NAME_PREFIX+amr_name+'_'+str(seq_length) in file_name:
			seq_file = file_name
	path_file = -1
	for file_name in path_info_file_names:
		if SEQ_NAME_PREFIX+amr_name+'_'+str(seq_length) in file_name:
			path_file = file_name
	return seq_file, path_file

def read_annotation_from_file(prokka_dir):
	"""
	to read annotations generated by Prokka and RGI for CAMI gold standard
	"""
	#Go over Prokka's output files and extract required information
	seq_info = []
	with open(prokka_dir+'mygenome.tsv', 'r') as tsvfile:
		reader = DictReader(tsvfile, delimiter='\t')
		for row in reader:
			mygene = row['gene'].strip()
			split_gene = mygene.split('_')
			if len(split_gene)==2 and split_gene[1].isnumeric():
				mygene = split_gene[0]
			gene_info = {'locus_tag':row['locus_tag'].strip(), 'gene':mygene,
						'length':row['length_bp'].strip(), 'contig': None,
						'product':row['product'].strip(),'start_pos':None, 'end_pos':None,
						'prokka_gene_name':mygene, 'RGI_prediction_type':None,
						'coverage':None, 'family': None, 'seq_value': None }
			seq_info.append(gene_info)
	counter = 0
	contig_name = ''
	with open(prokka_dir+'mygenome.tbl', 'r') as read_obj:
		for line in read_obj:
			if line.startswith('>'):
				contig_name = line.strip().split(' ')[1]
			elif line[0].isdigit():
				cells = line.split('\t')
				seq_info[counter]['start_pos'] = int(cells[0])
				seq_info[counter]['end_pos'] = int(cells[1])
				seq_info[counter]['contig'] = contig_name
				counter+=1
	#Read RGI output
	RGI_output_list = []
	if os.path.isfile(prokka_dir + 'rgi_output.txt'):
		with open(prokka_dir + 'rgi_output.txt', newline = '') as rgi_file:
			rgi_reader = csv.reader(rgi_file, delimiter='\t')
			next(rgi_reader)
			for row in rgi_reader:
				rgi_info = {'ORF_ID':row[0], 'gene':row[8].strip(),
				'prediction_type':row[5].strip(), 'best_identities':float(row[9]),
				'family':row[16].strip()}
				RGI_output_list.append(rgi_info)
	else:
		logging.error("ERROR: RGI didn't run successfully!")
		sys.exit()
	#incorporate RGI findings into Prokka's
	if RGI_output_list:
		for item in RGI_output_list:
			for gene_info in seq_info:
				if item['ORF_ID'].split(' ')[0]==gene_info['locus_tag']:
					gene_info['gene'] = item['gene']
					gene_info['RGI_prediction_type'] = item['prediction_type']
					gene_info['family'] = item['family']
					break

	return seq_info

def graph_extraction_main(params, gfa_file, graph_file, amr_seq_align_files):
	"""
	"""
	logging.info("Extracting neighborhood subgraphs ...")
	if not gfa_file:
		gfa_file = verify_file_existence(graph_file, params.gfa_file, \
				'please provide the address of the file containing the assembly graph')
	#remove paths from GFA file
	command = "sed -i '/^P/d' " + gfa_file
	os.system(command)
	if not os.path.exists(params.output_dir+SUBGRAPH_DIR_NAME):
		os.makedirs(params.output_dir+SUBGRAPH_DIR_NAME)
	if params.multi_processor:
		p_subgraph = partial(neighborhood_graph_extraction, gfa_file, params.graph_distance,
							params.output_dir+SUBGRAPH_DIR_NAME, params.BANDAGE_PATH,
							params.amr_identity_threshold, SEQ_NAME_PREFIX)
		with Pool(params.core_num) as p:
			p.map(p_subgraph, amr_seq_align_files)
	else:
		#for amr_file in amr_files:
		for amr_file in amr_seq_align_files:
			#output_name = 'ng_subgraph_'+os.path.splitext(os.path.basename(amr_file))[0]
			_ = neighborhood_graph_extraction(gfa_file, params.graph_distance,
				 				params.output_dir+SUBGRAPH_DIR_NAME, params.BANDAGE_PATH,
								params.amr_identity_threshold, SEQ_NAME_PREFIX, amr_file)


def seq_evaluation_main(params, amr_files, coverage_annotation_list,
				ref_up_info_lists, ref_down_info_lists, ref_amr_info_lists):
	"""
	"""
	logging.info("Neighborhood Evaluation ...")
	#storing summary metrics
	evaluation_dir = params.output_dir+EVAL_DIR+'/'+EVAL_DIR+'_'+str(params.seq_length)+'/'
	if not os.path.exists(evaluation_dir):
		try:
			os.makedirs(evaluation_dir)
		except OSError as exc:
			if exc.errno != errno.EEXIST:
				raise
			pass
	summary_file = evaluation_dir+'summaryMetrics_up_down_'+str(params.coverage_thr)+'_'+\
		datetime.datetime.now().strftime('%Y-%m-%d_%H-%M')+'.csv'
	with open(summary_file,'a') as fd:
		writer = csv.writer(fd)
		writer.writerow(['AMR', 'Unique_TP#', 'FP#', 'Unique_True#', 'found#','sensitivity', 'precision'])

	average_precision = 0
	average_sensitivity = 0
	for i, amr_file in enumerate(amr_files):
		amr_name = os.path.splitext(os.path.basename(amr_file))[0]
		# summary_file = params.output_dir+'evaluation/summaryMetrics_'+\
		# 	datetime.datetime.now().strftime('%Y-%m-%d_%H-%M')+'.csv'
		# with open(summary_file,'a') as fd:
		# 	writer = csv.writer(fd)
		# 	writer.writerow(['AMR', 'Unique_TP#', 'FP#', 'Unique_True#', 'found#','sensitivity', 'precision', 'Max Coverage unit diff for true cases', 'Min Coverage Unit diff for false cases'])
		if coverage_annotation_list and coverage_annotation_list[i]!="":
			sensitivity, precision = evaluate_sequences_up_down_based_on_coverage(
					amr_name, coverage_annotation_list[i], summary_file,
					ref_up_info_lists[i], ref_down_info_lists[i], ref_amr_info_lists[i])
			average_precision+=precision
			average_sensitivity+=sensitivity
			# sensitivity, precision = evaluate_sequences_based_on_coverage(output_name[1:], coverage_annotation,
			# 							evaluation_file, trimmed_annotation_info_name, amr_coverages)
			logging.info('For "'+amr_name+'": sensitivity= '+str(sensitivity)+' precision = '+ str(precision))
		elif not coverage_annotation_list:
			ref_up_info_list = []
			ref_down_info_list = []
			ref_amr_info_list = []
			annotate_dir = params.output_dir+ANNOTATION_DIR+'/'+ANNOTATION_DIR+'_'+str(params.seq_length)+'/annotation_'+amr_name+'_'+str(params.seq_length)
			coverage_annotation = annotate_dir+'/coverage_annotation_'+str(params.coverage_thr)+'_'+amr_name+'.csv'
			with open(coverage_annotation, 'r') as myfile:
				myreader = DictReader(myfile)
				old_seq = ''
				seq_info =[]
				for row in myreader:
					if not row['seq_name'].startswith('extracted'):
						gene_info = {'seq_name':row['seq_name'], 'seq_value':row['seq_value'],
						 			'gene':row['gene'], 'length':row['length'],
									'start_pos':int(row['start_pos']),'end_pos':int(row['end_pos']),
									'target_amr':row['target_amr']}
						cur_seq = row['seq_name']
						if cur_seq!=old_seq:
							if (seq_info):
								amr_found, up_info, down_info, amr_info = extract_up_down_from_csv_file(seq_info)
								if amr_found:
									ref_amr_info_list.append(amr_info)
									if up_info and not similar_seq_annotation_already_exist(up_info, ref_up_info_list):
										ref_up_info_list.append(up_info)
									if down_info and not similar_seq_annotation_already_exist(down_info, ref_down_info_list):
										ref_down_info_list.append(down_info)
							seq_info = []
							old_seq = cur_seq
						seq_info.append(gene_info)
				amr_found, up_info, down_info, amr_info = extract_up_down_from_csv_file(seq_info)
				if amr_found:
					ref_amr_info_list.append(amr_info)
					if up_info and not similar_seq_annotation_already_exist(up_info, ref_up_info_list):
						ref_up_info_list.append(up_info)
					if down_info and not similar_seq_annotation_already_exist(down_info, ref_down_info_list):
						ref_down_info_list.append(down_info)
			sensitivity, precision = evaluate_sequences_up_down_based_on_coverage(
					amr_name, coverage_annotation, summary_file,
					ref_up_info_list, ref_down_info_list, ref_amr_info_list)
			average_precision+=precision
			average_sensitivity+=sensitivity
			logging.info('For "'+amr_name+'": sensitivity= '+str(sensitivity)+' precision = '+ str(precision))

	return average_precision/len(amr_files), average_sensitivity/len(amr_files)

def seq_annotation_trim_main(params, amr_files, all_seq_info_lists, ref_seq_info_lists,
								ref_amr_info_lists, annotation_files, visualize = False):
	"""
	"""
	coverage_annotation_list = []
	for i, amr_file in enumerate(amr_files):
		amr_name = os.path.splitext(os.path.basename(amr_file))[0]
		#remove some extracted sequences based on coverage consistency
		annotate_dir = params.output_dir+ANNOTATION_DIR+'/'+ANNOTATION_DIR+'_'+str(params.seq_length)+'/annotation_'+amr_name+'_'+str(params.seq_length)
		coverage_annotation, remained_seqs = check_coverage_consistency_remove_rest_seq(\
								all_seq_info_lists[i], ref_seq_info_lists[i], params.output_dir,
								params.coverage_thr, amr_name, annotate_dir+'/')
		if visualize:
			# create an image presenting the annotations for all sequences
			if coverage_annotation!='':
				visual_annotation_csv = coverage_annotation
				visual_bound = remained_seqs + len(ref_amr_info_lists[i])
			else:
				visual_annotation_csv = annotation_files[i]
				visual_bound = len(all_seq_info_lists[i]) + len(ref_amr_info_lists[i])
			visual_annotation = annotate_dir+'/gene_comparison_'+str(params.coverage_thr)+'_'+amr_name+'.png'
			if visual_bound<20:
				visualize_annotation(visual_annotation_csv, output=visual_annotation)

		coverage_annotation_list.append(coverage_annotation)
	return coverage_annotation_list

def seq_annotation_main(params, seq_files, path_info_files, genome_amr_files,
						amr_files, ref_genome_files):
	"""
	"""
	logging.info("Neighborhood Annotation ...")
	if seq_files:
		neighborhood_files = seq_files
	else:
		neighborhood_files = extract_files(params.ng_seq_files, 'please provide the \
			address of the files containing all extracted sequences from AMR neighborhood \
			in the assembly graph')
	if path_info_files:
		nodes_info_files = path_info_files
	else:
		nodes_info_files = extract_files(params.ng_path_info_files, '')
	# contig_file = verify_file_existence(contigs_file, params.contig_file,
	# 		'please provide the address of the file containing contigs after assembly')
	if params.artificial_amr_insertion:
		if genome_amr_files:
			genome_files = genome_amr_files
		elif not params.genome_amr_files:
			logging.error('ERROR: please provide the address of the files containing genome after AMR insertion')
			#print('ERROR: please provide the address of the files containing genome after AMR insertion')
			sys.exit()
		else:
			genome_files = params.genome_amr_files
	else:
		genome_files = ref_genome_files

	#If it's for a CAMI sample, we need to do some processing on ref annotation beforehand
	if params.ref_CAMI_genome_available and CAMI_REF_NEIGHBORHOOD_METHOD in [1, 2]:
		amr_family = extract_amr_family_info(AMR_FAMILY_INFO)
		ref_seq_info = read_annotation_from_file(params.output_dir+'prokka_dir/')

	ref_up_info_lists = []
	ref_down_info_lists = []
	ref_amr_info_lists = []
	all_seq_info_lists = []
	ref_seq_info_lists = []
	annotation_files = []
	for amr_file in amr_files:
		amr_name = os.path.splitext(os.path.basename(amr_file))[0]
		cami_info = None
		if params.ref_CAMI_genome_available and CAMI_REF_NEIGHBORHOOD_METHOD in [1, 2]:
			#find the family of this gene
			myfamily = ''
			for k, gene_list in enumerate(amr_family['gene_list']):
				if amr_name.lower() in gene_list:
					myfamily = amr_family['family'][k]
					break
			if myfamily=='':
				logging.error('ERROR: no family was found for '+amr_name)
				break
			logging.debug('amr_family: '+myfamily)
			cami_info = (myfamily, ref_seq_info, params.ref_CAMI_file)
		elif params.ref_CAMI_genome_available and CAMI_REF_NEIGHBORHOOD_METHOD==3:
			cami_info = ('', [], params.ref_CAMI_file)
			#Create DB from genome files
			command = 'makeblastdb -in '+params.ref_CAMI_file +' -parse_seqids -dbtype nucl'
			os.system(command)

		neighborhood_file, nodes_info_file = find_corrsponding_seq_path_file(amr_name, neighborhood_files,
												nodes_info_files, params.seq_length)
		if neighborhood_file == -1:
			logging.error('no sequence file for the corresponding amr file was found!')
			#print('no sequence file for the corresponding amr file was found!')
			sys.exit()
		if params.multi_processor:
			all_seq_info_list, ref_seq_info_list, ref_up_info_list, ref_down_info_list,\
			ref_amr_info_list, annotation_file=\
				neighborhood_annotation_parallel(amr_file, genome_files,neighborhood_file,
					nodes_info_file, params.seq_length,
					params.output_dir, params.PROKKA_COMMAND_PREFIX,params.use_RGI,
					params.RGI_include_loose, '_'+amr_name,
					params.amr_identity_threshold, params.ref_genomes_available,
					params.core_num, cami_info)
		else:
			all_seq_info_list, ref_seq_info_list, ref_up_info_list, ref_down_info_list,\
			ref_amr_info_list, annotation_file =\
				neighborhood_annotation(amr_file, genome_files,neighborhood_file,
					nodes_info_file, params.seq_length,
					params.output_dir, params.PROKKA_COMMAND_PREFIX,params.use_RGI,
					params.RGI_include_loose, '_'+amr_name,
					params.amr_identity_threshold, params.ref_genomes_available,
					cami_info)

		ref_up_info_lists.append(ref_up_info_list)
		ref_down_info_lists.append(ref_down_info_list)
		ref_amr_info_lists.append(ref_amr_info_list)
		all_seq_info_lists.append(all_seq_info_list)
		ref_seq_info_lists.append(ref_seq_info_list)
		annotation_files.append(annotation_file)

	return all_seq_info_lists, ref_seq_info_lists, ref_up_info_lists,\
			ref_down_info_lists, ref_amr_info_lists, annotation_files

def sequence_neighborhood_main(params, gfa_file, graph_file, amr_seq_align_files):
	"""
	"""
	seq_files = []
	path_info_files = []
	logging.info("Extracting neighborhood sequences with length = %s", params.seq_length)
	if not gfa_file:
		gfa_file = verify_file_existence(graph_file, params.gfa_file, \
				'please provide the address of the file containing the assembly graph')
	#remove paths from GFA file
	command = "sed -i '/^P/d' " + gfa_file
	os.system(command)
	sequence_dir = params.output_dir+SEQ_DIR_NAME+'/'+SEQ_DIR_NAME+'_'+str(params.seq_length)+'/'
	if not os.path.exists(sequence_dir):
		os.makedirs(sequence_dir)
	if params.multi_processor:
		p_extraction = partial(neighborhood_sequence_extraction, gfa_file, params.seq_length,
							sequence_dir, params.BANDAGE_PATH,
							params.amr_identity_threshold, SEQ_NAME_PREFIX,
							params.path_node_threshold , params.path_seq_len_percent_threshod,
							params.max_kmer_size, params.assembler)
		with Pool(params.core_num) as p:
		#with Pool(4) as p:
			lists = p.map(p_extraction, amr_seq_align_files)
			#lists = p.map(p_extraction, amr_files)
		seq_files, path_info_files = zip(*lists)
	else:
		#for amr_file in amr_files:
		for amr_file in amr_seq_align_files:
			#output_name = SEQ_NAME_PREFIX+os.path.splitext(os.path.basename(amr_file))[0]
			seq_file, path_info_file = neighborhood_sequence_extraction(gfa_file, params.seq_length,
								sequence_dir, params.BANDAGE_PATH,
								params.amr_identity_threshold, SEQ_NAME_PREFIX,
								params.path_node_threshold , params.path_seq_len_percent_threshod,
								params.max_kmer_size, params.assembler, amr_file)
			if seq_file:
				path_info_files.append(path_info_file)
				seq_files.append(seq_file)

	return seq_files, path_info_files

def main(params):
	logging.info("Startting the pipeline ...")
	#Validate task values
	task_list = validate_task_values(params.task)
	if params.artificial_amr_insertion and params.find_amr_genes:
		logging.error("variables 'artificial_amr_insertion' and 'find_amr_genes' cannot be True at the same run!")
		#print("variables 'artificial_amr_insertion' and 'find_amr_genes' cannot be True at the same run!")
		sys.exit()

	graph_file =""
	# read1 = ""
	# read2 = ""
	metagenome_file = ""
	genome_amr_files = []
	path_info_files = []
	seq_files = []
	reads = ""
	gfa_file = None
	amr_files = []
	amr_align_files = []
	ref_genome_files = []
	if params.ref_genomes_available:
		ref_genome_files = extract_files(params.ref_genome_files, 'please provide the address of genome files')
		#Create DB from genome files
		for genome_file in ref_genome_files:
			command = 'makeblastdb -in '+genome_file +' -parse_seqids -dbtype nucl'
			os.system(command)


	if not params.find_amr_genes:
		amr_files = extract_files(params.amr_files, 'please provide the address of the AMR gene(s)')
		align_dir = os.path.abspath(os.path.join(os.path.dirname(params.amr_files),'..',AMR_ALIGN_DIR))
		align_files = extract_files(align_dir, '')
		#remove align files which are in the overlap list
		if align_files:
			for amr_file in amr_files:
				found_it= False
				amr_name = os.path.splitext(os.path.basename(amr_file))[0]
				for align_file in align_files:
					if os.path.basename(align_file).startswith(amr_name+'_align'):
						found_it=True
						amr_align_files.append(align_file)
						break
				if not found_it:
					import pdb; pdb.set_trace()
	if params.ref_genomes_available and params.artificial_amr_insertion:
		if not amr_files:
			logging.error("ERROR: no AMR file has been provided!")
			#print("ERROR: no AMR file has been provided!")
			sys.exit()
		# for now, we only accept one AMR gene to be inserted in the genomes
		if len(amr_files)>1:
			logging.error("ERROR: for the artificial AMR insertion, we can't accept more than one AMR gene sequence!")
			#print("ERROR: for the artificial AMR insertion, we can't accept more than one AMR gene sequence!")
			sys.exit()

	if params.ref_genomes_available and Pipeline_tasks.metagenome_creation.value in task_list:
		logging.info("Creating the metagenome sample ...")
		if params.artificial_amr_insertion:
			genome_amr_files, metagenome_file = create_metagenome_with_amr_insertion(ref_genome_files,
						params.number_of_insertions, params.insertion_type, params.insertion_locations,
						amr_files[0], params.output_dir)
		else:
			metagenome_file = create_metagenome_file(ref_genome_files, params.output_dir)

	if Pipeline_tasks.read_simulation.value in task_list:
		logging.info("Simulating reads ...")
		meta_file = verify_file_existence(metagenome_file, params.metagenome_file, \
					'please provide the address of metagenome file')
		reads = simulate_reads(meta_file, params.read_length, params.ART_PATH)

	if Pipeline_tasks.assembly.value in task_list:
		logging.info("Assembly ...")
		# read1_file = verify_file_existence(read1, params.read1, \
		# 	'please provide the address of both paired end reads')
		# read2_file = verify_file_existence(read2, params.read2, \
		# 	'please provide the address of both paired end reads')
		#import pdb; pdb.set_trace()
		read_files = verify_file_existence(reads, params.reads, \
			'please provide the address of paired end reads')
		spade_output = params.output_dir + params.assembler_output_dir
		graph_file, contigs_file = do_assembly(read_files, params.SPADES_PATH,
										spade_output, params.spades_thread_num,
										params.spades_error_correction)
		# graph_file, contigs_file = do_assembly(read1_file, read2_file, params.SPADES_PATH,
		# 								spade_output, params.spades_thread_num,
		# 								params.spades_error_correction)

	if params.find_amr_genes:
		logging.info("Finding AMR genes in the assembly graph ...")
		gfa_file = verify_file_existence(graph_file, params.gfa_file, \
				'please provide the address of the file containing the assembly graph')
		if params.multi_processor:
			amr_files, amr_align_files = find_all_amr_in_graph_parallel(gfa_file, params.output_dir,
							params.PATH_PREFIX+ALL_AMR_SEQUENCES, params.BANDAGE_PATH,
							params.amr_identity_threshold, params.core_num)
		else:
			amr_files, amr_align_files = find_all_amr_in_graph(gfa_file, params.output_dir,
							params.PATH_PREFIX+ALL_AMR_SEQUENCES, params.BANDAGE_PATH,
							params.amr_identity_threshold)
	if not amr_files:
		logging.info("No AMR gene was found!")
		sys.exit()
	send_amr_align_files = False
	if amr_align_files and len(amr_align_files)==len(amr_files):
		send_amr_align_files = True
	#create pairs of seq and align files
	amr_seq_align_files = []
	for i, amr_file in enumerate(amr_files):
		if send_amr_align_files:
			amr_seq_align_files.append((amr_file, amr_align_files[i]))
		else:
			amr_seq_align_files.append((amr_file, ''))

	# if Pipeline_tasks.graph_neighborhood.value in task_list:
	# 	graph_extraction_main(params, gfa_file, graph_file, amr_seq_align_files)

	# if Pipeline_tasks.sequence_neighborhood.value in task_list:

	if Pipeline_tasks.sequence_neighborhood.value in task_list:
		if MULTIPLE_SEQ_LENGTH:
			for seq_len in [500, 1000, 1500, 2000, 2500, 3000, 3500, 4000, 4500, 5000]:
				params.seq_length = seq_len
				seq_files, path_info_files = sequence_neighborhood_main(params, gfa_file,
						graph_file, amr_seq_align_files)
		else:
			seq_files, path_info_files = sequence_neighborhood_main(params, gfa_file,
					graph_file, amr_seq_align_files)

	coverage_annotation_list = []
	ref_up_info_lists = []
	ref_down_info_lists = []
	ref_amr_info_lists = []
	if Pipeline_tasks.neighborhood_annotation.value in task_list:
		all_seq_info_lists, ref_seq_info_lists, ref_up_info_lists,\
		ref_down_info_lists, ref_amr_info_lists, annotation_file_list =\
			seq_annotation_main(params, seq_files, path_info_files,
					genome_amr_files, amr_files, ref_genome_files)
		if MULTIPLE_COV_THR:
			evaluation_dir = params.output_dir+EVAL_DIR+'/'+EVAL_DIR+'_'+str(params.seq_length)+'/'
			if not os.path.exists(evaluation_dir):
				try:
					os.makedirs(evaluation_dir)
				except OSError as exc:
					if exc.errno != errno.EEXIST:
						raise
					pass
			coverage_evaluation_file = evaluation_dir+'coverage_evaluation_'+\
				datetime.datetime.now().strftime('%Y-%m-%d_%H-%M')+'.csv'
			with open(coverage_evaluation_file,'a') as fd:
				writer = csv.writer(fd)
				writer.writerow(['coverage_thr', 'precision', 'sensitivity'])
				data = []
				for cov_thr in range(1, 50):
					params.coverage_thr = cov_thr
					coverage_annotation_list = seq_annotation_trim_main(params, amr_files,\
					 		all_seq_info_lists, ref_seq_info_lists, ref_amr_info_lists, annotation_file_list, False)
					if Pipeline_tasks.neighborhood_evaluation.value in task_list:
						precision, sensitivity = seq_evaluation_main(params, amr_files, coverage_annotation_list,
									ref_up_info_lists, ref_down_info_lists, ref_amr_info_lists)
						data.append([cov_thr, precision, 'Precision'])
						data.append([cov_thr, sensitivity, 'Sensitivity'])
						writer.writerow([cov_thr, precision, sensitivity])
			df = pd.DataFrame(data, columns = ['cov_thr', 'value', 'type'])
			sns.scatterplot(data=df, x = 'cov_thr', y = 'value', hue='type', style='type')
			plt.show()
		else:
			coverage_annotation_list = seq_annotation_trim_main(params, amr_files,\
				all_seq_info_lists, ref_seq_info_lists, ref_amr_info_lists, annotation_file_list, True)

	if not MULTIPLE_COV_THR and Pipeline_tasks.neighborhood_evaluation.value in task_list:
		seq_evaluation_main(params, amr_files, coverage_annotation_list,
					ref_up_info_lists, ref_down_info_lists, ref_amr_info_lists)

	logging.info("All Done!")

if __name__=="__main__":

	import params

	text = 'This code is used to find the context of a given AMR gene'
	parser = argparse.ArgumentParser(description=text)
	parser.add_argument('--task', nargs="+", default=params.task,
		help="which task would you like to do?\
		For the entire pipeline choose "+str(Pipeline_tasks.all.value)+"; otherwise\
		either provide a number representing one of the following tasks or two numbers\
		to denote the start and end tasks (and of course all tasks in the middle will be run).\n \
		Here is the list:\nmetagenome_creation = "+str(Pipeline_tasks.metagenome_creation.value)+\
		"\nread_simulation = "+str(Pipeline_tasks.read_simulation.value)+
		"\nassembly = "+str(Pipeline_tasks.assembly.value)+
		"\ngraph_neighborhood = "+str(Pipeline_tasks.graph_neighborhood.value)+\
		"\nsequence_neighborhood = "+str(Pipeline_tasks.sequence_neighborhood.value)+\
		"\nneighborhood_annotation = "+str(Pipeline_tasks.neighborhood_annotation.value)+\
		"\nneighborhood_evaluation = "+str(Pipeline_tasks.neighborhood_evaluation.value))
	parser.add_argument('--amr_files','-A', type=str, default = params.amr_files,
		help = 'the path of the file(s) containing the AMR gene sequence(s)')
	parser.add_argument('--ref_genome_files', nargs="+", default=params.ref_genome_files,
		help = 'the ddress of reference genomes that AMR genome will be inserted in them')
	parser.add_argument('--ref_CAMI_file', nargs="+", default=params.ref_CAMI_file,
		help = 'the ddress of reference genomes that AMR genome will be inserted in them for CAMI')
	parser.add_argument('--output_dir', '-O', type = str, default=params.output_dir,
		help = 'the output dir to store metagenome file and reads')
	parser.add_argument('--number_of_insertions', type = int, default=params.number_of_insertions,
		help = 'the number of genomes generated by inserting AMR in different locations of reference genome')
	parser.add_argument('--insertion_type', type =Insertion_type , default=params.insertion_type,
		help = 'Should insertion locations be selected randomly (1) or using some defined values (2)')
	parser.add_argument('--insertion_locations', nargs="+", default = params.insertion_locations,
		help = 'list of predefined insertion locations to insert AMR in reference genome\
			if you chose multiple reference genomes first type the list of insertion locations\
			for the first genome followed by the list for the second and ...\
			please make sure the number of locations for each genome matches the\
			corresponding one specified in --num_insertion argument')
	parser.add_argument('--metagenome_file', type = str, default= params.metagenome_file,
		help = 'the address of metagenome file')
	parser.add_argument('--read_length',type=int, default=params.read_length,
		help = 'the length of simulated reads can be either 150 or 250')
	parser.add_argument('--spades_thread_num',type=int, default=params.spades_thread_num,
		help = 'the number of threads used for MetaSPAdes')
	parser.add_argument('--assembler_output_dir',type=str, default=params.assembler_output_dir,
		help = 'the output dir to store MetaSPAdes results')
	parser.add_argument('--graph_distance', '-D', type = int, default=params.graph_distance,
		help = 'the maximum distance of neighborhood nodes to be extracted from the AMR gene')
	parser.add_argument('--seq_length', '-L', type = int, default=params.seq_length,
		help = 'the length of AMR gene\'s neighbourhood to be extracted')
	parser.add_argument('--ng_seq_files', nargs="+", default = params.ng_seq_files,
		help = 'the address of the files containing all extracted neighborhood sequences in assembly graph')
	parser.add_argument('--ng_path_info_files', nargs="+", default = params.ng_path_info_files,
		help = 'the address of the files containing all path information for extracted neighborhood sequences in assembly graph')
	parser.add_argument('--gfa_file', type = str, default = params.gfa_file,
		help = 'the address of the file for assembly graph')
	parser.add_argument('--contig_file', type = str, default = params.contig_file,
		help = 'the address of the file containing contigs after assembly')
	parser.add_argument('--genome_amr_files', nargs="+", default = params.genome_amr_files,
		help = 'the address of the files containing genome after AMR insertion')
	# parser.add_argument('--read1', type = str, default = params.read1,
	# 	help = 'the address of the files containing the first file in paired-end reads')
	# parser.add_argument('--read2', type = str, default = params.read2,
	# 	help = 'the address of the files containing the second file in paired-end reads')
	parser.add_argument('--reads', type = check_reads, default = params.reads,
		help = 'the address of the files containing paired-end reads')
	parser.add_argument('--spades_error_correction', type = str2bool, default = params.spades_error_correction,
		help = 'Whether to turn on or off error correction in MetaSPAdes')
	parser.add_argument('--use_RGI', type = str2bool, default = params.use_RGI,
		help = 'Whether to contribute RGI annotation in Prokka result')
	parser.add_argument('--RGI_include_loose', type = str2bool, default = params.RGI_include_loose,
		help = 'Whether to include loose cases in RGI result')
	parser.add_argument('--find_amr_genes', type = str2bool, default = params.find_amr_genes,
		help = 'Whether to assume the AMR genes (in metagenome) are known or to look for them in assembly graph')
	parser.add_argument('--artificial_amr_insertion', type = str2bool, default = params.artificial_amr_insertion,
		help = 'Whether to insert the AMR gene in genomes artificially')
	parser.add_argument('--amr_identity_threshold', type = int, default = params.amr_identity_threshold,
		help = 'the threshold used for amr alignment: a hit is returned if identity/coverage >= threshold')
	parser.add_argument('--path_node_threshold', type = int, default = params.path_node_threshold,
		help = 'the threshold used for recursive pre_path and post_path search as long as the length of the path is less that this threshold')
	parser.add_argument('--path_seq_len_percent_threshod', type = int, default = params.path_seq_len_percent_threshod,
		help = 'the threshold used for recursive pre_seq and post_seq until we have this percentage of the required length\
		 after which we just extract from the longest neighbor')
	parser.add_argument('--ref_genomes_available', type = str2bool, default = params.ref_genomes_available,
		help = 'Whether we have access to reference genome(s)')
	parser.add_argument('--ref_CAMI_genome_available', type = str2bool, default = params.ref_CAMI_genome_available,
		help = 'Whether we have access to reference genome(s) for CAMI')
	parser.add_argument('--multi_processor', type = str2bool, default = params.multi_processor,
		help = 'Whether to use multi processors for parallel programming')
	parser.add_argument('--core_num', type = int, default = params.core_num,
		help = 'the number of cores used in case of parallel programming')
	parser.add_argument('--coverage_thr', type = int, default = params.coverage_thr,
		help = 'coverage threshold to check if an annotated gene is truely AMR neighbor or just a falose positive')

	args = parser.parse_args()

	#updating params values
	#parameters = {}
	params.task = args.task
	params.amr_files = args.amr_files
	params.ref_genome_files = args.ref_genome_files
	params.ref_CAMI_file = args.ref_CAMI_file
	params.output_dir = args.output_dir
	params.number_of_insertions = args.number_of_insertions
	params.insertion_type = args.insertion_type
	params.insertion_locations = args.insertion_locations
	params.read_length = args.read_length
	params.spades_thread_num = args.spades_thread_num
	params.assembler_output_dir = args.assembler_output_dir
	params.graph_distance = args.graph_distance
	params.seq_length = args.seq_length
	params.ng_seq_files = args.ng_seq_files
	params.ng_path_info_files =  args.ng_path_info_files
	params.gfa_file = args.gfa_file
	params.contig_file = args.contig_file
	params.genome_amr_files = args.genome_amr_files
	# params.read1 = args.read1
	# params.read2 = args.read2
	params.reads = args.reads
	params.spades_error_correction = args.spades_error_correction
	params.use_RGI = args.use_RGI
	params.RGI_include_loose = args.RGI_include_loose
	params.find_amr_genes = args.find_amr_genes
	params.metagenome_file =  args.metagenome_file
	params.artificial_amr_insertion = args.artificial_amr_insertion
	params.amr_identity_threshold = args.amr_identity_threshold
	params.path_node_threshold = args.path_node_threshold
	params.path_seq_len_percent_threshod = args.path_seq_len_percent_threshod
	params.ref_genomes_available = args.ref_genomes_available
	params.ref_CAMI_genome_available = args.ref_CAMI_genome_available
	params.coverage_thr = args.coverage_thr
	params.multi_processor = args.multi_processor
	params.core_num = args.core_num

	log_name = 'logger_'+datetime.datetime.now().strftime('%Y-%m-%d_%H-%M')+'.log'
	initialize_logger(params.output_dir, log_name)
	logging.info(str(params.__dict__))
	#print_params(params)

	main(params)
