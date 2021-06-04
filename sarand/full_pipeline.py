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
import errno
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
import math

from params import Pipeline_tasks, Insertion_type, Assembler_name
from extract_neighborhood import neighborhood_graph_extraction, neighborhood_sequence_extraction,\
							extract_amr_neighborhood_in_ref_genome
from find_seq_in_contigs import find_sequence_match
from annotation_visualization import visualize_annotation
from utils import initialize_logger, check_reads, str2bool, print_params, verify_file_existence,\
			retrieve_AMR, extract_files, create_fasta_file, extract_amr_names_from_alignment_files,\
			annotate_sequence, split_up_down_info, unnamed_genes_are_siginificantly_similar,\
			seqs_annotation_are_identical, similar_seq_annotation_already_exist,\
			read_ref_annotations_from_db, extract_up_down_from_csv_file, amr_name_from_comment,\
			amr_name_from_title, retreive_original_amr_name, extract_name_from_file_name,\
			restricted_amr_name_from_modified_name, extract_info_from_overlap_file,\
			read_info_from_overlap_ref_files, extract_unique_align_files,\
			read_path_info_from_align_file, concatenate_files, read_path_info_from_align_file_with_multiple_amrs,\
			extract_path_info_for_amrs, compare_two_sequences, split_up_down_seq,\
			delete_lines_started_with
from find_amrs_in_sample import find_annotate_amrs_in_ref

ASSEMBLY_FILE = 'assembly_graph_with_scaffolds.gfa'
CONTIG_FILE = 'contigs.fasta'
#ALL_AMR_SEQUENCES ='nucleotide_fasta_protein_homolog_model_without_efflux_without_space.fasta'
AMR_DIR_NAME = 'AMR_info/'
AMR_SEQ_DIR = 'sequences/'
AMR_ALIGN_DIR = 'alignments/'
AMR_OVERLAP_FILE = 'overlaps.txt'
SUBGRAPH_DIR_NAME = 'subgraphs/'
SEQ_DIR_NAME = 'sequences_info'
SEQ_NAME_PREFIX = 'ng_sequences_'
ANNOTATION_DIR = 'annotations'
EVAL_DIR = 'evaluation'
NOT_FOUND_FILE = 'not_found_amrs_in_graph.txt'

#To run the code for a list of sequence neighborhood length rather than just one length
#the default seq length is 1000
MULTIPLE_SEQ_LENGTH = False
#To Do the sequence evaluation rather than annotation evaluation
SEQ_EVAL = False

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
		import pdb; pdb.set_trace()
		sys.exit()

	valid_task_values = [item.value for item in Pipeline_tasks]
	for task in tasks:
		if int(task) not in valid_task_values:
			logging.error("ERROR: invalid task number(s)!\n" + task_error_message)
			import pdb; pdb.set_trace()
			sys.exit()

	if len(tasks)==2 and int(tasks[0])>int(tasks[1]):
		logging.error("ERROR: The first task number should be smaller than the second task\
		 in the list!\n" + task_error_message)
		import pdb; pdb.set_trace()
		sys.exit()

	if len(tasks)==1 and int(tasks[0])==Pipeline_tasks.all.value:
		return valid_task_values
	if len(tasks)==1:
		return [int(tasks[0])]
	for task in list(range(int(tasks[0]), int(tasks[1])+1)):
		task_list.append(task)
	return task_list

def find_amrs_not_in_graph(ref_amr_files, amr_names):
	"""
	compare the list of ref amr files with the amr names available in the graph
	to find the ones not available in the graph
	Parameter:
		ref_amr_files: the list of amr files available in ref genomes
		amr_names: the list of amr names available in the graph
	Return:
		the list of amr names not found in the graph
	"""
	not_found = []
	if len(ref_amr_files)==len(amr_names):
		return not_found
	elif len(ref_amr_files)<len(amr_names):
		logging.error("We are supposed to only work on the AMRs found on ref genomes and nothing beyond that")
		import pdb; pdb.set_trace()
	ref_amr_names = [extract_name_from_file_name(e) for e in ref_amr_files]
	for amr_name in amr_names:
		if amr_name not in ref_amr_names:
			not_found.append(amr_name)

	return not_found

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
		use_RGI: if True, RGI has been used to annotate AMRs
		found: if True, the annotation info has already found in other annotated sequences
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
		import pdb; pdb.set_trace()
		sys.exit()
	if params.insertion_type == Insertion_type.assigned and \
		len(params.insertion_locations)!= sum(params.number_of_insertions):
		logging.error("ERROR: Please specify the location of insertion for all insertions OR \
			if you prefer them to be chosen randomely choose 1 for --insertion_type")
		import pdb; pdb.set_trace()
		sys.exit()

	if not os.path.exists(output_dir):
		os.makedirs(output_dir)
	if (os.path.isfile(output_dir+metagenome_file)):
		logging.error('ERROR: A file named ' + metagenome_file + ' already exits in '+output_dir)
		import pdb; pdb.set_trace()
		sys.exit()

	insert_counter = 0
	genome_files = []
	amr_seq, _ = retrieve_AMR(amr_file)
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
		import pdb; pdb.set_trace()
		sys.exit()
	logging.info("Running ART: "+command)
	os.system(command)

	read1 = metagenome_file_name + '_1.fq'
	read2 = metagenome_file_name + '_2.fq'
	return [read1, read2]

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
		command = spades_path + ' -1 '+reads[0]+' -2 '+reads[1]+' --meta '
	else:
		command = spades_path + ' --12 '+reads + ' --meta '
	if not error_correction:
		command+='--only-assembler '
	command += '--threads ' + str(thread_num) + ' -o '+ output_dir
	logging.info("Running MetaSPAdes: "+command)
	os.system(command)

	#remove paths from GFA file
	gfa_file = output_dir + '/'+ ASSEMBLY_FILE
	delete_lines_started_with('P', gfa_file)

	contig_file = output_dir + '/' + CONTIG_FILE

	return gfa_file, contig_file

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
			found = True
			break

	return found

def read_path_info_file(path_info_file):
	"""
	To read the csv file containing the list of nodes and their coverage representing
	in an extracted sequence
	Parameter:
		path_info_file: the csv file
	Return:
		stored info in a list
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
	Calculate the coverage of genes available in a sequence based on the coverage
	of the nodes representing them
	Parameter:
		seq_info_list: the list of genes (and their info) annotated in a sequence
		path_info:	the list of nodes (and their info) representing the sequence
	Return:
		the list of calculated gene-coverages
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
			logging.error("ERROR: no nodes were found for this gene!!!")
			import pdb; pdb.set_trace()
			sys.exit()
	return coverage_list

def this_is_target_amr_row(amr_name, gene, heads, member_lists):
	"""
	To find out if the current gene is the target AMR gene or not!
	We first compare gene's name with the name of the amr;
	However, since sometimes the gene representing amr is annotated with a different name,
	we also check all members of the group represented by our AMR
	Parameter:
		amr_name: the name of target AMR
		gene: the name of the gene, we are currently processing
		heads: the list of heads of the groups
		members_lists: the lists of members of each group
	Return:
		True if gene is representing AMR and False otherwise
	"""
	gene_name = amr_name_from_title(gene)
	if amr_name.strip()==gene_name.strip():
		return True
	for i, item in enumerate(heads):
		if item==amr_name:
			for member in member_lists[i]:
				if member.strip()==gene.strip():
					return True
	return False

def find_target_amr_in_seqvalue_and_return_coverage(seq_info):
	"""
	To find the target AMR in a sequence and its covearge:
	We first look at the extracted sequence as the part of the sequence in lower-case
	represents the AMR; then, we look at the list of annotated genes and find the
	gene that has the most overlap with the lower-case area and finally find its coverage
	Parameters:
		seq_info: the details of an extracted sequences and its annotation
	Return:
		the coverage of the AMR, its index in seq_info and False if we can find a
		gene in the lower-case area of the sequence; otherwise it returns True
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

def check_coverage_consistency_remove_rest_seq(seq_info_list_input,
								coverage_thr, amr_name, annotate_dir):
	"""
	To compare the coverage of each gene in the neighborhood with that of the AMR
	and remove the gene (and any gene before that in case of upstream and any gene
	after that in case of downstream) if the difference between this gene's coverage
	and ame coverage is more than covearge_thr
	Parameters:
		seq_info_list_input: the list of all sequences and their info including annotations
		coverage_thr:	the threshold used to compare the gene covearge and amr covearge
		amr_name:	the name of target AMR
		annotate_dir:	the directory in which annotation info are stored
	Return:
	"""
	seq_info_list = copy.deepcopy(seq_info_list_input)
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
			if gene_info['target_amr']=='yes':
				amr_coverages.append(coverage)
				amr_indeces.append(gene_counter)
				found_amr = True
				break
		if not found_amr:
			amr_coverage, amr_index, error = find_target_amr_in_seqvalue_and_return_coverage(seq_info)
			if error:
				logging.info("ERROR: no target amr was found for "+ str(seq_info)+" regarding "+amr_name)
				import pdb; pdb.set_trace()
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
		#write extracted sequences with consistent coverage
		for seq_info in remained_seqs:
			for gene_info in seq_info:
				writer.writerow([gene_info['seq_name'], gene_info['seq_value'],
							len(gene_info['seq_value']),
							gene_info['gene'], gene_info['coverage'],
							gene_info['length'], gene_info['start_pos'],
							gene_info['end_pos'], gene_info['target_amr']])

	return coverage_annotation, len(remained_seqs)

def exists_in_gene_list(gene_info_lists, gene_info):
	"""
	To check if a gene (and its annotation info) are available in a list of
	annoated sequences. we compare the name of the genes.
	In case that the target gene has no name, we call another
	function to check if its sequence is significantly similar to another unnamed
	gene in the list.
	Parameters:
		gene_info_lists:	the list of annotated sequences
		gene_info:	the gene we are searching for
	Return:
		True if it can find the gene in the list; False otherwise
	"""
	for gene_info_list in gene_info_lists:
		for item in gene_info_list:
			if gene_info['gene']!='' and item['gene']==gene_info['gene']:
				return True
			elif gene_info['gene']=='' and item['gene']==gene_info['gene']:
				return unnamed_genes_are_siginificantly_similar(item, gene_info, threshold = 90)
	return False

def evaluate_sequences_up_down_based_on_coverage(amr_name, coverage_annotation, summary_file,
												ref_up_info_list, ref_down_info_list,
												ref_amr_info_list, assembler):
	"""
	To compare the annotated sequences extracted from the graph with the ones from
	the ref genomes. For any AMR, each graph extracted upstream (downstream)
	annotation will be compared with the ref extracted upstream (downstream),
	and if all genes are the same in both annoations it will be considered a true positive!
	Parameters:
		amr_name: target AMR
		coverage_annotation: the file containing annotation of all extracted sequences for this AMR
		summary_file: the file to store the results (precision and sensitivity for AMRs)
		ref_up_info_list: the list of annotation of ref extracted up streams
		ref_down_info_list: the list of annotation of ref extracted down streams
		ref_amr_info_list: the list of annotation of amr in ref
		assembler: the name of assembler used to generated the assembly graph
	Return:
		calculated precision and sensitivity
	"""
	ref_len =  len(ref_up_info_list)+len(ref_down_info_list)
	original_amr_name = retreive_original_amr_name(amr_name)
	#if amr was not found in the ref genomes
	if ref_len==0 and len(ref_amr_info_list)==0:
	#This shouldn't happen when ref_genome is available!!!
		with open(summary_file,'a') as fd:
			writer = csv.writer(fd)
			writer.writerow([original_amr_name, 0, 0, 0, len(up_info_list)+len(down_info_list),-1, -1])
		logging.error(amr_name+" was not found in the ref genomes!!!")
		import pdb; pdb.set_trace()
		return -1, -1
	#If AMR was not found in the contig list
	if coverage_annotation=="":
		with open(summary_file,'a') as fd:
			writer = csv.writer(fd)
			writer.writerow([original_amr_name, 0, 0, ref_len, 0,0, 0])
		return 0, 0
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
				if row['coverage']=='' and assembler!=Assembler_name.metacherchant:
					logging.error("ERROR: Coverage information are not available for "+amr_name)
					import pdb; pdb.set_trace()
					sys.exit()
				#else:
				gene_info = {'seq_name':row['seq_name'], 'seq_value':row['seq_value'],
				 			'gene':row['gene'], 'length':row['length'],
							'start_pos':int(row['start_pos']),'end_pos':int(row['end_pos']),
							'target_amr':row['target_amr']}
				cur_seq = row['seq_name']
				if cur_seq!=old_seq:
					if (seq_info):
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
		amr_found, up_info, down_info, amr_info = extract_up_down_from_csv_file(seq_info)
		if amr_found:
			amr_info_list.append(amr_info)
			if up_info and not similar_seq_annotation_already_exist(up_info, up_info_list):
				up_info_list.append(up_info)
			if down_info and not similar_seq_annotation_already_exist(down_info, down_info_list):
				down_info_list.append(down_info)

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
	if found_cases_len == 0 and ref_len == 0:
		precision = 1
	else:
		precision = 0 if found_cases_len==0 else round(1 - float(false_positive)/found_cases_len, 2)

	with open(summary_file,'a') as fd:
		writer = csv.writer(fd)
		writer.writerow([original_amr_name, unique_tp, false_positive, ref_len, found_cases_len,
						sensitivity, precision])

	return sensitivity, precision

def evaluate_sequences_gene_based_on_coverage(amr_name, coverage_annotation, summary_file,
											annotation_file, amr_coverages):
	"""
	To compare the annotated sequences extracted from the graph with the ones from
	the ref genomes. For any AMR, each graph extracted gene will be compared with
	the corresponding ref extracted gene,and if they are the same in both annoations
	it will be considered a true positive!
	Parameters:
		amr_name: target AMR
		coverage_annotation: the file containing annotation of all extracted sequences for this AMR
		summary_file: the file to store the results (precision and sensitivity for AMRs)
		annotation_file: the file containting ref info
		amr_coverage: the coverage of AMR
	Return:
		calculated precision and sensitivity

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
					import pdb; pdb.set_trace()
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
		logging.error("ERROR: inconsistency between length of amr_coverages and seq_info_list2")
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
		original_amr_name = retreive_original_amr_name(amr_name)
		writer.writerow([original_amr_name, true_positive, false_positive, len(ref_info_list), total,
						sensitivity, precision, max_true, min_false if found_false else ''])

	return sensitivity, precision

def extract_seq_annotation(annotate_dir, prokka_prefix, use_RGI, RGI_include_loose,
								seq_pair):
	"""
	The function used in parallel anntation to call the function for annotating a sequence
	Parameters:
		annotate_dir: the directory to store annotation output
		prokka_prefix: the first part of prokka command probably is nonempty only when docker is used
		use_RGI: if True we want to call RGI for annotating AMRs
		RGI_include_loose: if True loose mode is used
		seq_pair: the index and the value of a sequence to be annotated
	Return:
		the list of annotated genes
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
	Parameters:
		sequence: the sequence that has been annotated
		seq_info: the annotation info
		up_info_list: the annotation of the upstream sequences in the reference genomes
		down_info_list: the annotation of downstream sequences in the reference genomes
	Return:
		modified annotation info for the sequence
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
		else:
			# added by 1 because in string indecesstarts from 0
			diff = max((amr_start+1-start), 0)+max((end - (amr_end+1)), 0)
			if ((1-(float(diff)/(end-start)))*100)>overlap_thr:
				found = True
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
	return seq_info

def add_coverage_to_info(seq_info, up_info, down_info, amr_info, coverage_list):
	"""
	To add gene coverage to the gene information dictionary
	Parameters:
		seq_info: annotation info for a sequence including the information for annotated genes
		up_info: the annotation info for the upstream
		down_info: the annotation info for the downstream
		amr_info: the annotation info for the target AMR
		coverage_list: the list of coverage calculated for genes in seq_info
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
	amr_seq, _ = retrieve_AMR('411.fasta')
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
			#genome_name = os.path.splitext(os.path.basename(genome_file))[0]
			genome_name = extract_name_from_file_name(genome_file)
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
	To extract the neighborhood of an AMR from ref genomes and annotate it
	Parameters:
		amr_file: the file containing AMR sequence
		output_name: the name of AMR used for naming files
		seq_length: the neighborhood sequence length
		ref_genome_files: the list of files containing ref genomes
		amr_threshold: the threshold for identity and coverage to align AMR using blastn
		annotate_dir: the output dir to store annotation results
		prokka_prefix: the prefix for prokka command, mostly used for docker
		use_RGI: if True RGI is used to annotate AMRs
		RGI_include_loose: if True include loose mode of RGI
		annotation_writer: file to add annotation detail
		trimmed_annotation_writer: another file to add annotation detail
		gene_file: file to store annotation info in terms of the name of genes
		product_file: file to store annotation info in terms of the name of products
	Return:
		annotation info for upstream, downsytream, AMR and all
	"""
	ref_up_info_list = []
	ref_down_info_list = []
	ref_amr_info_list = []
	#not necessarily unique values are stored here
	ref_seq_info_list = []
	amr_seq, _ = retrieve_AMR(amr_file)
	ref_file_name = annotate_dir+'/ref_neighborhood_sequences'+output_name+'_'+\
		str(seq_length)+'_'+datetime.datetime.now().strftime('%Y-%m-%d_%H-%M')+'.txt'
	ref_file = open(ref_file_name, 'w')
	for genome_file in ref_genome_files:
		seq_list, _ = extract_amr_neighborhood_in_ref_genome(amr_seq, genome_file, seq_length, amr_threshold)
		#AMR sequence might be in multiple places in ref genome so we have a list of
		#extracted neighborhood sequences instead of one sequence
		for index, seq in enumerate(seq_list):
			#genome_name = os.path.splitext(os.path.basename(genome_file))[0]
			genome_name = extract_name_from_file_name(genome_file)
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
	To read a part of a contig from a fasta file
	Parameters:
		contig_name: the name of the contig
		up_limit: the starting point of the target sequence
		down_limit: the end point of the target sequence
		ref_CAMI_file: the file in which the search is done
	Return:
		the extracted sequence or Error if it can't be found
	"""
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
		#gene = gene_info['gene'].strip().replace(' ','_').replace("'",';').replace('/', ']')
		#gene_name = ''.join(e for e in gene if e.isalpha() or e.isnumeric() or e=='_' or e=='-')
		gene_name = amr_name_from_title(gene_info['gene'])
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
	To check if the AMR can be found the list annotated genes of the same contig after index position
	Parameters:
		amr_name: the name of AMR
		ref_seq_info: the annotation list
		index: the position of interest in  ref_seq_info
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
	amr_seq, _ = retrieve_AMR(amr_file)
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

def extract_graph_seqs_annotation_parallel(amr_name, path_info_file, neighborhood_seq_file,
					annotate_dir, core_num, ref_genomes_available, prokka_prefix,
					use_RGI, RGI_include_loose,ref_up_info_list, ref_down_info_list,
					annotation_writer, trimmed_annotation_writer, gene_file,
					product_file, error_file):
	"""
	To annotate neighborhood sequences of AMR extracted from the graph in parallel
	Parameters:
		amr_name: the name of AMR
		path_info_file: the information of nodes representing the extracted sequence
		neighborhood_seq_file: the file containing extracted sequences
		annotate_dir: the directory to sore annotation results
		core_num: the number of core for parallel processing
		ref_genomes_available: if True reference genomes are available
		prokka_prefix: the prefix in prokka command mostly used in docker
		use_RGI: if True RGI is used for AMR annotation
		RGI_include_loose: if True use loose mode in RGI
		ref_up_info_list: the annotation detail for upstream extracted from ref genomes
		ref_down_info_list: the annotation detail for downstream extracted from ref genomes
		annotation_writer: the file to store annotation results
		trimmed_annotation_writer: the file to store unique annotation results
		gene_file: the file to store gene nams in annotation
		product_file: the file to store product name in annotation
		error_file: the file to store errors
	Return:
		the list of annotated genes and their details
	"""
	error_writer = open(error_file, 'a')
	#Read path_info from the file
	path_info_list = []
	if path_info_file!=-1:
		path_info_list = read_path_info_file(path_info_file)
	#find the list of all extracted sequences
	logging.info('Reading '+ neighborhood_seq_file + ' for '+ amr_name)
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
			error_writer.write(amr_name+' annotation not found! '+" seq_info: "+str(seq_info)+'\n')
			continue
			#import pdb; pdb.set_trace()
		#calculate coverage for the genes available in the annotation
		coverage_list = []
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
			coverage = coverage_list[j] if coverage_list else -1
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
	if not all_seq_info_list:
		error_writer.write(amr_name+' no annotation was found in the graph.\n')
	error_writer.close()
	return all_seq_info_list

def neighborhood_annotation_parallel(amr_name, neighborhood_seq_file,
								path_info_file, seq_length,
								ref_up_info_list, ref_down_info_list,
								output_dir, prokka_prefix, use_RGI = True,
								RGI_include_loose = False, output_name ='',
								amr_threshold = 95, ref_genomes_available = True,
								core_num = 4):
	"""
	To annotate reference genomes (a piece extracted around the AMR gene) as well as
		extracted neighborhood sequences from assembly graph, summarize the results
		in a couple of formats and visualize them
	Parameters:
		amr_name:	the name of target AMR
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
	logging.debug('Started annotation for '+amr_name)
	# initializing required files and directories
	annotate_dir = output_dir+ANNOTATION_DIR+'/'+ANNOTATION_DIR+'_'+str(seq_length)+'/annotation'+output_name+'_'+str(seq_length)
	if os.path.exists(annotate_dir):
		try:
			shutil.rmtree(annotate_dir)
		except OSError as e:
			logging.error("Error: %s - %s." % (e.filename, e.strerror))
	os.makedirs(annotate_dir)
	error_file =output_dir+ANNOTATION_DIR+'/'+ANNOTATION_DIR+'_'+str(seq_length) + "/not_found_annotation_amrs_in_graph.txt"
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

	#annotate the sequences extraced from assembly graph
	all_seq_info_list = extract_graph_seqs_annotation_parallel(amr_name, path_info_file, neighborhood_seq_file,
									annotate_dir, core_num, ref_genomes_available,
									prokka_prefix, use_RGI, RGI_include_loose,
									ref_up_info_list, ref_down_info_list,
									annotation_writer, trimmed_annotation_writer,
									gene_file, product_file, error_file)
	logging.info("NOTE: The comparison of neighborhood sequences are available in " +\
	 		annotation_detail_name+", "+gene_file_name+", "+product_file_name)
	annotation_detail.close()
	trimmed_annotation_info.close()
	gene_file.close()
	product_file.close()

	return all_seq_info_list, trimmed_annotation_info_name

def extract_graph_seqs_annotation(amr_name, path_info_file, neighborhood_seq_file,
					annotate_dir, ref_genomes_available, prokka_prefix,
					use_RGI, RGI_include_loose,ref_up_info_list, ref_down_info_list,
					annotation_writer, trimmed_annotation_writer, gene_file, product_file):
	"""
	To annotate neighborhood sequences of AMR extracted from the graph
	Parameters:
		amr_name: the name of AMR
		path_info_file: the information of nodes representing the extracted sequence
		neighborhood_seq_file: the file containing extracted sequences
		annotate_dir: the directory to sore annotation results
		ref_genomes_available: if True reference genomes are available
		prokka_prefix: the prefix in prokka command mostly used in docker
		use_RGI: if True RGI is used for AMR annotation
		RGI_include_loose: if True use loose mode in RGI
		ref_up_info_list: the annotation detail for upstream extracted from ref genomes
		ref_down_info_list: the annotation detail for downstream extracted from ref genomes
		annotation_writer: the file to store annotation results
		trimmed_annotation_writer: the file to store unique annotation results
		gene_file: the file to store gene nams in annotation
		product_file: the file to store product name in annotation
		error_file: the file to store errors
	Return:
		the list of annotated genes and their details
	"""
	error_file ="not_found_annotation_amrs_in_graph.txt"
	error_writer = open(error_file, 'a')
	counter = 1
	all_seq_info_list =[]
	#Read path_info from the file
	path_info_list = []
	if path_info_file!=-1 and path_info_file!='':
		path_info_list = read_path_info_file(path_info_file)
	logging.info('Reading '+ neighborhood_seq_file + ' for '+ amr_name)
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
				error_writer.write(amr_name+' annotation not found! '+" seq_info: "+str(seq_info)+'\n')
				continue
				#import pdb; pdb.set_trace()
			#calculate the coverage of annotated genes
			coverage_list = []
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
				coverage = coverage_list[j] if coverage_list else -1
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
	if not all_seq_info_list:
		error_writer.write(amr_name+' no annotation was found in the graph.\n')
	error_writer.close()
	return all_seq_info_list

def neighborhood_annotation(amr_name, neighborhood_seq_file,
								path_info_file, seq_length,
								ref_up_info_list, ref_down_info_list,
								output_dir, prokka_prefix, use_RGI = True,
								RGI_include_loose = False, output_name ='',
								amr_threshold = 95, ref_genomes_available = True):
	"""
	To annotate reference genomes (a piece extracted around the AMR gene) as well as
		extracted neighborhood sequences from assembly graph, summarize the results
		in a couple of formats and visualize them.
	Parameters:
		amr_name:	the name of target AMR
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
	logging.debug('Started annotation for '+amr_name)
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
	#annotate the sequences extraced from assembly graph
	all_seq_info_list = extract_graph_seqs_annotation(amr_name, path_info_file, neighborhood_seq_file,
					annotate_dir, ref_genomes_available, prokka_prefix,
					use_RGI, RGI_include_loose,ref_up_info_list, ref_down_info_list,
					annotation_writer, trimmed_annotation_writer, gene_file, product_file)
	annotation_detail.close()
	trimmed_annotation_info.close()
	gene_file.close()
	product_file.close()
	logging.info("NOTE: The comparison of neighborhood sequences are available in " +\
	 		annotation_detail_name+", "+gene_file_name+", "+product_file_name)

	return all_seq_info_list, trimmed_annotation_info_name

def is_there_amr_in_graph_check_longer_paths_incrementally(amr_name, gfa_file, output_dir,
						bandage_path, threshold, amr_file, MAX_PATH_NODES = 50):
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
	#amr_name = os.path.splitext(os.path.basename(amr_file))[0]
	amr_name = extract_name_from_file_name(amr_file)
	logging.info('Checking if AMR "'+amr_name+'" exists in the assembly graph...')
	output_name=output_dir+amr_name+'_align_'+datetime.datetime.now().strftime('%Y-%m-%d_%H-%M-%S')
	found = False
	path_nodes = 10
	while (not found) and (path_nodes < MAX_PATH_NODES):
		if os.path.isfile(output_name+'.tsv'):
			os.remove(output_name+'.tsv')
		command = bandage_path +' querypaths ' + gfa_file+' '+amr_file+' '+output_name +\
		' --pathnodes '+str(path_nodes)
		os.system(command)
		found, paths_info = read_path_info_from_align_file(output_name+".tsv", threshold)
		path_nodes+=1

	if not found:
		os.remove(output_name+'.tsv')
		logging.debug(amr_name+' not found in the graph!')
	else:
		logging.debug(amr_name+' found!!!')
	return found, paths_info, output_name+".tsv"

def is_there_amr_in_graph(gfa_file, output_dir, bandage_path, threshold, amr_file):
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
	amr_name = extract_name_from_file_name(amr_file)
	logging.info('Checking if AMR "'+amr_name+'" exists in the assembly graph...')
	output_name=output_dir+amr_name+'_align_'+datetime.datetime.now().strftime('%Y-%m-%d_%H-%M-%S')
	if os.path.isfile(output_name+'.tsv'):
		os.remove(output_name+'.tsv')
	command = bandage_path +' querypaths '+gfa_file+' '+amr_file+' '+output_name + ' --pathnodes 50'
	os.system(command)

	found, paths_info = read_path_info_from_align_file(output_name+".tsv", threshold)
	if not found:
		os.remove(output_name+'.tsv')
		logging.debug(amr_name+' not found in the graph!')
	else:
		logging.debug(amr_name+' found!')
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
						if i not in id_list:
							id_list.append(i)
						break
			if found:
				break
	if len(id_list)==len(new_paths):
		return True, id_list
	return False, None

def are_there_amrs_in_graph(gfa_file, output_dir, bandage_path, threshold, amr_object):
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
	cat_file, amr_files = amr_object
	amr_names = [extract_name_from_file_name(e) for e in amr_files]
	logging.info('Checking if AMRs "'+str(amr_names)+'" exists in the assembly graph...')
	output_name=output_dir+extract_name_from_file_name(cat_file)+'_align_'+datetime.datetime.now().strftime('%Y-%m-%d_%H-%M-%S')
	if os.path.isfile(output_name+'.tsv'):
		os.remove(output_name+'.tsv')
	command = bandage_path +' querypaths '+gfa_file+' '+cat_file+' '+output_name + ' --pathnodes 50'
	os.system(command)

	paths_info_list = read_path_info_from_align_file_with_multiple_amrs(output_name+".tsv", threshold)

	return paths_info_list

def find_amrs_alignment_in_graph_parallel(gfa_file, output_dir, amr_files, bandage_path,
											amr_threshold, core_num):
	"""
	To go over a list of AMR sequences and run bandage+blast
	to check if any of them exists in the assembly graph (gfa_file) and
	return their alignment files if they are found in the graph
	Parameters:
		gfa_file: the address of the assembly graph
		output_dir: the address of the output directory
		amr_files: the list of address of AMR files
		bandage_path: the address of bandage executable file
		amr_threshold: the threshold for coverage and identity
		core_num: the number of used cores
	Return:

	"""
	align_dir = output_dir+AMR_DIR_NAME+AMR_ALIGN_DIR
	if not os.path.exists(align_dir):
		os.makedirs(align_dir)
	else:
		try:
			shutil.rmtree(align_dir)
		except OSError as e:
			logging.error("Error: %s - %s." % (e.filename, e.strerror))
		os.makedirs(align_dir)
	#generate the groups and store the group of each amr
	group_num = 2
	amr_group_id = collections.defaultdict(list)
	amr_groups_num = min(group_num*core_num, len(amr_files))
	amr_file_groups = [[] for i in range(amr_groups_num)]
	for i, amr_file in enumerate(amr_files):
		id = i % (amr_groups_num)
		amr_file_groups[id].append(amr_file)
		amr_group_id[amr_file] = id
	#concatenate files of each group into a single file
	amr_objects = []
	for i, file_group in enumerate(amr_file_groups):
		cat_file = output_dir+AMR_DIR_NAME+ 'amr_group_'+str(i)+'.fasta'
		concatenate_files(file_group, cat_file)
		amr_objects.append((cat_file, file_group))
	#find AMRs parallel
	p_find_amr_align = partial(are_there_amrs_in_graph, gfa_file, align_dir,
							bandage_path, amr_threshold)
	with Pool(core_num) as p:
		paths_info_group_list = p.map(p_find_amr_align, amr_objects)

	#write the list of AMRs not found
	found_writer = open(output_dir+AMR_DIR_NAME+NOT_FOUND_FILE, 'w')
	not_found_amr_names = []
	unique_amr_files = []
	unique_amr_infos = []
	unique_amr_paths = []
	align_files = []
	for i, amr_file in enumerate(amr_files):
		id = amr_group_id[amr_file]
		restricted_amr_name = extract_name_from_file_name(amr_file)
		seq, amr_name = retrieve_AMR(amr_file)
		if restricted_amr_name in paths_info_group_list[id]:
			logging.debug(amr_name+' was found!')
			path_info = paths_info_group_list[id][restricted_amr_name]
			overlap, amr_ids =  amr_path_overlap(unique_amr_paths, path_info, len(seq[:-1]))
			if not overlap:
				unique_amr_files.append(amr_file)
				amr_info = {'name':amr_name, 'overlap_list':[]}
				unique_amr_infos.append(amr_info)
				unique_amr_paths.append(path_info)
			else:
				if len(amr_ids)>1:
					logging.error("an AMR has overlap with more than one group")
					import pdb; pdb.set_trace()
				# add this AMR to the right group of AMRs all having overlaps
				for id in amr_ids:
					if amr_name not in unique_amr_infos[id]['overlap_list']:
						unique_amr_infos[id]['overlap_list'].append(amr_name)
		else:
			not_found_amr_names.append(amr_name)
			found_writer.write(amr_name+'\n')
	# write the list of groups in each all AMRs have overlaped paths into a file
	overlap_file_name = output_dir+AMR_DIR_NAME+AMR_OVERLAP_FILE
	overlap_file = open(overlap_file_name, 'w')
	for amr_info in unique_amr_infos:
		overlap_file.write(amr_info['name']+":")
		if amr_info['overlap_list']:
			overlap_file.write(', '.join(e for e in amr_info['overlap_list']))
			overlap_file.write("\n")
		else:
			overlap_file.write("\n")
	overlap_file.close()
	found_writer.close()
	return unique_amr_files, not_found_amr_names, unique_amr_paths

def find_amrs_alignment_in_graph(gfa_file, output_dir, amr_files, bandage_path, amr_threshold):
	"""
	To go over a list of AMR sequences and run bandage+blast
	to check if any of them exists in the assembly graph (gfa_file) and
	return their alignment files if they are found in the graph
	Parameters:
		gfa_file: the address of the assembly graph
		output_dir: the address of the output directory
		amr_files: the list of address of AMR files
		bandage_path: the address of bandage executable file
		amr_threshold: the threshold for coverage and identity
	Return:

	"""
	align_dir = output_dir+AMR_DIR_NAME+AMR_ALIGN_DIR
	if not os.path.exists(align_dir):
		os.makedirs(align_dir)
	else:
		try:
			shutil.rmtree(align_dir)
		except OSError as e:
			logging.error("Error: %s - %s." % (e.filename, e.strerror))
		os.makedirs(align_dir)

	#write the list of AMRs not found
	found_writer = open(output_dir+AMR_DIR_NAME+NOT_FOUND_FILE, 'w')
	not_found_amr_names = []
	unique_amr_files = []
	unique_amr_infos = []
	unique_amr_paths = []
	for amr_file in amr_files:
		#amr_name = os.path.splitext(os.path.basename(amr_file))[0]
		seq, amr_name = retrieve_AMR(amr_file)
		found, paths_info, tsv_file = is_there_amr_in_graph( gfa_file, align_dir,
										bandage_path, amr_threshold, amr_file)
		if found:
			seq = seq[:-1]
			overlap, amr_ids =  amr_path_overlap(unique_amr_paths, paths_info, len(seq))
			if not overlap:
				unique_amr_files.append(amr_file)
				amr_info = {'name':amr_name, 'overlap_list':[]}
				unique_amr_infos.append(amr_info)
				unique_amr_paths.append(paths_info)
			else:
				# add this AMR to the right group of AMRs all having overlaps
				if len(amr_ids)>1:
					logging.error("an AMR has overlap with more than one group")
					import pdb; pdb.set_trace()
				for id in amr_ids:
					if amr_name not in unique_amr_infos[id]['overlap_list']:
						unique_amr_infos[id]['overlap_list'].append(amr_name)
		else:
			not_found_amr_names.append(amr_name)
			found_writer.write(amr_name+'\n')
	# write the list of groups in each all AMRs have overlaped paths into a file
	overlap_file_name = output_dir+AMR_DIR_NAME+AMR_OVERLAP_FILE
	overlap_file = open(overlap_file_name, 'w')
	for amr_info in unique_amr_infos:
		overlap_file.write(amr_info['name']+":")
		if amr_info['overlap_list']:
			#overlap_file.write(amr_info['name']+":\n")
			overlap_file.write(', '.join(e for e in amr_info['overlap_list']))
			overlap_file.write("\n")
		else:
			overlap_file.write("\n")
	overlap_file.close()
	found_writer.close()
	return unique_amr_files, not_ound_amr_names, unique_amr_paths

def process_amr_group_and_find(gfa_file, align_dir, output_dir, bandage_path,
										amr_threshold, amr_object):
	"""
	Read a group of AMRs, write them into a single file and call bandage+blast for it
	to identify them in the graph
	Parameters:
		gfa_file: the file containing the assembly graph
		align_dir: the directory for storing alignment info
		output_dir: the directory to store the list of AMRs in a single file
		bandage_path: the path for bandage
		amr_threshor: the threshold for identity and coverage used in alignment
		amr_object: the list of AMRs and their ids
	Return:
		the alignment info for AMRs
	"""
	g_id, amr_group = amr_object
	#read info of the group into a single file
	cat_file = output_dir+AMR_DIR_NAME+ 'amr_group_'+str(g_id)+'.fasta'
	file_group = []
	with open(cat_file, 'w') as writer:
		for amr_info in amr_group:
			amr_seq, amr_title = amr_info
			writer.write(amr_title)
			writer.write(amr_seq)
			amr_name1 = amr_name_from_comment(amr_title)
			amr_file_name = restricted_amr_name_from_modified_name(amr_name1)
			file_group.append(amr_file_name+'.fasta')

	#Run Bandage+BLAST
	p_find_amr_align = are_there_amrs_in_graph(gfa_file, align_dir,
										bandage_path, amr_threshold, (cat_file, file_group))
	#Remove temporary AMR file
	if os.path.isfile(cat_file):
		os.remove(cat_file)
	return p_find_amr_align

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
		core_num: the number of used cores
	Return:

	"""
	align_dir = output_dir+AMR_DIR_NAME+AMR_ALIGN_DIR
	if not os.path.exists(align_dir):
		os.makedirs(align_dir)

	#generate the groups and store the group of each amr
	group_num = 5
	amr_group_id = collections.defaultdict(list)
	amr_file_groups = [[] for i in range(group_num*core_num)]
	amr_title = ''
	amr_seq_title_list = []
	#Read AMR sequences one by one
	amr_counter = 0
	with open(amr_sequences_file) as fp:
		for line in fp:
			if line.startswith('>'):
				amr_title = line
				continue
			amr_name = amr_name_from_comment(amr_title[:-1])
			amr_seq_title_list.append((line, amr_title))
			id = amr_counter % (group_num*core_num)
			amr_file_groups[id].append((line, amr_title))
			amr_group_id[amr_name] = id
			amr_counter+=1

	amr_objects = [(i, e) for i,e in enumerate(amr_file_groups)]
	#parallel run Bandage+BLAST
	p_find_amr = partial(process_amr_group_and_find, gfa_file, align_dir,
							output_dir, bandage_path, amr_threshold)
	with Pool(core_num) as p:
		paths_info_group_list = p.map(p_find_amr, amr_objects)


	unique_amr_seqs = []
	unique_amr_infos = []
	unique_amr_paths = []
	#process the result of parallel processes
	for i, amr_object in enumerate(amr_seq_title_list):
		amr_name = amr_name_from_comment(amr_object[1])
		id = amr_group_id[amr_name]
		restricted_amr_name = restricted_amr_name_from_modified_name(amr_name)
		if restricted_amr_name in paths_info_group_list[id]:
			logging.debug(amr_name + ' was found!')
			path_info = paths_info_group_list[id][restricted_amr_name]
			overlap, amr_ids =  amr_path_overlap(unique_amr_paths, path_info,
													len(amr_object[0])-1)
			if not overlap:
				unique_amr_seqs.append(amr_object[0])
				amr_info = {'name':amr_object[1], 'overlap_list':[]}
				unique_amr_infos.append(amr_info)
				unique_amr_paths.append(path_info)
			else:
				if len(amr_ids)>1:
					logging.error("an AMR has overlap with more than one group")
					import pdb; pdb.set_trace()
				# add this AMR to the right group of AMRs all having overlaps
				for id in amr_ids:
					if amr_name not in unique_amr_infos[id]['overlap_list']:
						unique_amr_infos[id]['overlap_list'].append(amr_name)

	# write information (the sequence of found AMRs that don't have overlaped paths with others
	# + the list of groups in each all AMRs have overlaped paths) into files
	AMR_dir = output_dir+AMR_DIR_NAME+AMR_SEQ_DIR
	if not os.path.exists(AMR_dir):
		os.makedirs(AMR_dir)
	overlap_file_name = output_dir+AMR_DIR_NAME+AMR_OVERLAP_FILE
	overlap_file = open(overlap_file_name, 'w')
	unique_amr_files = []
	for i, seq in enumerate(unique_amr_seqs):
		amr_name = amr_name_from_comment(unique_amr_infos[i]['name'])
		restricted_amr_name = restricted_amr_name_from_modified_name(amr_name)
		amr_file = create_fasta_file(seq, AMR_dir, unique_amr_infos[i]['name'], restricted_amr_name)
		unique_amr_files.append(amr_file)
		overlap_file.write(amr_name+":")
		if unique_amr_infos[i]['overlap_list']:
			overlap_file.write(', '.join(e for e in unique_amr_infos[i]['overlap_list']))
			overlap_file.write("\n")
		else:
			overlap_file.write("\n")
	overlap_file.close()

	return unique_amr_files, unique_amr_paths

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
	#found_amr_names = []
	unique_amr_seqs = []
	unique_amr_infos = []
	unique_amr_paths = []
	#Read AMR sequences one by one
	with open(amr_sequences_file) as fp:
		for line in fp:
			if line.startswith('>'):
				amr_title = line
				continue
			#create a fasta file for it
			amr_file = create_fasta_file(line, output_dir)
			#Run Bandage+BLAST
			amr_name = amr_name_from_comment(amr_title)
			found, paths_info, tsv_file = is_there_amr_in_graph(gfa_file, align_dir,
											bandage_path, amr_threshold, amr_file)
			if found:
				#found_amr_names.append(amr_name)
				overlap, amr_ids =  amr_path_overlap(unique_amr_paths, paths_info,
										len(line)-1)
				if not overlap:
					unique_amr_seqs.append(line)
					amr_info = {'name':amr_title, 'overlap_list':[]}
					unique_amr_infos.append(amr_info)
					unique_amr_paths.append(paths_info)
				else:
					if len(amr_ids)>1:
						logging.error("an AMR has overlap with more than one group")
						import pdb; pdb.set_trace()
					# add this AMR to the right group of AMRs all having overlaps
					for id in amr_ids:
						if amr_name not in unique_amr_infos[id]['overlap_list']:
							unique_amr_infos[id]['overlap_list'].append(amr_name)
	# write information (the sequence of found AMRs that don't have overlaped paths with others
	# + the list of groups in each all AMRs have overlaped paths) into files
	AMR_dir = output_dir+AMR_DIR_NAME+AMR_SEQ_DIR
	if not os.path.exists(AMR_dir):
		os.makedirs(AMR_dir)
	overlap_file_name = output_dir+AMR_DIR_NAME+AMR_OVERLAP_FILE
	overlap_file = open(overlap_file_name, 'w')
	unique_amr_files = []
	for i, seq in enumerate(unique_amr_seqs):
		amr_name = amr_name_from_comment(unique_amr_infos[i]['name'])
		restricted_amr_name = restricted_amr_name_from_modified_name(amr_name)
		amr_file = create_fasta_file(seq, AMR_dir, unique_amr_infos[i]['name'], restricted_amr_name)
		unique_amr_files.append(amr_file)
		overlap_file.write(amr_name+":")
		if unique_amr_infos[i]['overlap_list']:
			#overlap_file.write(amr_name+":\n")
			overlap_file.write(', '.join(e for e in unique_amr_infos[i]['overlap_list']))
			overlap_file.write("\n")
		else:
			overlap_file.write("\n")
	overlap_file.close()

	return unique_amr_files, unique_amr_paths

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
		import pdb; pdb.set_trace()
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

def evaluate_pure_sequences(amr_name, neighborhood_seq_file, ref_up_seq_list,
							ref_down_seq_list, evaluation_csv, summary_file,threshold):
	"""
	to compare the sequences extracted from the graph from those of the ref genomes directly
	rather than comparing their annotated version.
	upstream and downstreams are compared separately
	Parameters:
		amr_name: the name of AMR
		neighborhood_seq_file: the file containing extracted sequences from the graph
		ref_up_seq_list: extracted upstreams from ref genomes
		ref_down_seq_list: extracted downstreams from ref genomes
		evaluation_csv: the file to store the details of comparision
		summary_file: the file to store the preciosn and sensitivity for each AMR
		threshold: the threshold used for identity and coverage to consider two sequences similar enough
	Return:
		precision and sensitivity
	"""
	output_dir = os.path.dirname(evaluation_csv)+'/'
	ref_len =  len(ref_up_seq_list)+len(ref_down_seq_list)
	#read sequences from neighborhood_seq_file
	up_seq_list = []
	down_seq_list = []
	if neighborhood_seq_file!='':
		with open(neighborhood_seq_file, 'r') as read_obj:
			for line in read_obj:
				if line.startswith('>') or line.startswith('Path') or line.startswith('The'):
					continue
				up_seq, down_seq = split_up_down_seq(line[:-1])
				if up_seq!='' and up_seq not in up_seq_list:
					up_seq_list.append(up_seq)
				if down_seq!='' and down_seq not in down_seq_list:
					down_seq_list.append(down_seq)
	#if no seq neighborhood was available in the ref genomes
	if ref_len==0:
		with open(summary_file,'a') as fd:
			writer = csv.writer(fd)
			writer.writerow([amr_name, 0, 0, 0, len(up_seq_list)+len(down_seq_list),-1, -1])
		logging.error(amr_name+" was not found in the ref genomes!!!")
		import pdb; pdb.set_trace()
		return -1, -1
	#If AMR was not found in the contig list
	if neighborhood_seq_file=='':
		with open(summary_file,'a') as fd:
			writer = csv.writer(fd)
			writer.writerow([amr_name, 0, 0, ref_len, 0,0, 0])
		return 0, 0
	#store identity and coverage info
	#import pdb; pdb.set_trace()
	for up_seq in up_seq_list:
		for ref_up_seq in ref_up_seq_list:
			blast_file_name = compare_two_sequences(ref_up_seq, up_seq, output_dir,
					threshold, switch_allowed = False, return_file = True)
			any_row_found = False
			with open(blast_file_name, 'r') as file1:
				myfile = csv.reader(file1)
				identity = ''
				coverage = ''
				for row in myfile:
					any_row_found = True
					iden = float(row[2])
					cov = float(row[3])/len(ref_up_seq)*100
					if iden>=20 and cov>=20:
						identity = identity + str(iden)+','
						coverage = coverage + str(cov)+','
			#if not any_row_found:
			#	import pdb; pdb.set_trace()
			with open(evaluation_csv, 'a') as ed:
				eval_writer = csv.writer(ed)
				eval_writer.writerow([amr_name, 'up_stream', up_seq, ref_up_seq, identity[:-1], coverage[:-1]])
	for down_seq in down_seq_list:
		for ref_down_seq in ref_down_seq_list:
			blast_file_name = compare_two_sequences(ref_down_seq, down_seq, output_dir,
					threshold, switch_allowed = False, return_file = True)
			any_row_found = False
			with open(blast_file_name, 'r') as file1:
				myfile = csv.reader(file1)
				identity = ''
				coverage = ''
				for row in myfile:
					any_row_found = True
					iden = float(row[2])
					cov = float(row[3])/len(ref_down_seq)*100
					if iden>=20 and cov>=20:
						identity = identity + str(iden)+','
						coverage = coverage + str(cov)+','
			#if not any_row_found:
			#	import pdb; pdb.set_trace()
			with open(evaluation_csv, 'a') as ed:
				eval_writer = csv.writer(ed)
				eval_writer.writerow([amr_name, 'down_stream', down_seq, ref_down_seq, identity[:-1], coverage[:-1]])

	#find the number of unique true-positives, all false positives, total found cases, all unique true cases
	unique_tp = 0
	for ref_seq in ref_up_seq_list:
		for seq in up_seq_list:
			if compare_two_sequences(ref_seq, seq, output_dir, threshold, switch_allowed = False):
				unique_tp+=1
				break
	for ref_seq in ref_down_seq_list:
		for seq in down_seq_list:
			if compare_two_sequences(ref_seq, seq, output_dir, threshold, switch_allowed = False):
				unique_tp+=1
				break
	sensitivity = 1 if ref_len==0 else round(float(unique_tp)/ref_len, 2)

	false_positive = 0
	for seq in up_seq_list:
		found_similar_seq = False
		for ref_seq in ref_up_seq_list:
			if compare_two_sequences(ref_seq, seq, output_dir, threshold, switch_allowed = False):
				found_similar_seq = True
				break
		if not found_similar_seq:
			false_positive+=1
	for seq in down_seq_list:
		found_similar_seq = False
		for ref_seq in ref_down_seq_list:
			if compare_two_sequences(ref_seq, seq, output_dir, threshold, switch_allowed = False):
				found_similar_seq = True
				break
		if not found_similar_seq:
			false_positive+=1
	found_cases_len = len(up_seq_list)+len(down_seq_list)
	if found_cases_len == 0 and ref_len == 0:
		precision = 1
	else:
		precision = 0 if found_cases_len==0 else round(1 - float(false_positive)/found_cases_len, 2)

	with open(summary_file,'a') as fd:
		writer = csv.writer(fd)
		writer.writerow([amr_name, unique_tp, false_positive, ref_len, found_cases_len,
						sensitivity, precision])


	return precision, sensitivity

def graph_extraction_main(params, gfa_file, graph_file, amr_seq_align_files):
	"""
	The core function to extract the neighborhood graph for AMRs
	Parameters:
		params: list of parameters set in params.py
		gfa_file: the file containing the assembly graph
		graph_file: the file containing the assembly graph generated in the program
		amr_seq_align_files: the list of AMRs and their alignment info
	"""
	logging.info("Extracting neighborhood subgraphs ...")
	if not gfa_file:
		gfa_file = verify_file_existence(graph_file, params.gfa_file, \
				'please provide the address of the file containing the assembly graph')
	#remove paths from GFA file
	delete_lines_started_with('P', gfa_file)
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


def annotation_evaluation_main(params, amr_files, coverage_annotation_list,
				ref_up_info_lists, ref_down_info_lists, ref_amr_info_lists, not_found_amr_names):
	"""
	the core function to evaluate annotations of all AMRs
	Parameters:
		params: the list of parameetrs imported from params.py
		amr_files: the list of files containing AMRs
		coverage_annotation_list: the list of csv files containing annotation info
		ref_up_info_lists: the list of annotations of extracted upstream sequences from ref genomes
		ref_down_info_lists: the list of annotations of extracted downstream sequences from ref genomes
		ref_amr_info_lists: the list of annotations of extracted amr sequences from ref genomes
		not_found_amr_names: the list of AMRs not found in the graph which are available in the ref genome
	Return:
		average precision and average sensitivity
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
	#extract the list of AMR groups
	overlap_file_name = params.output_dir+AMR_DIR_NAME+AMR_OVERLAP_FILE
	heads, member_lists, _ = extract_info_from_overlap_file(overlap_file_name)
	#to find the annotation of ref genomes for all AMRs
	df = pd.read_csv(params.ref_ng_annotations_file, skipinitialspace=True,  keep_default_na=False)
	amr_groups = df.groupby('target_amr')
	#Find the summary metrics for all AMRs
	average_precision = 0
	average_sensitivity = 0
	#Go over the list of unique AMRs (the heads in overlap file + some other not in overlap list)
	for i, amr_file in enumerate(amr_files):
		restricted_amr_name = extract_name_from_file_name(amr_file)
		_, amr_name = retrieve_AMR(amr_file)
		if coverage_annotation_list and coverage_annotation_list[i]!="":
			sensitivity, precision = evaluate_sequences_up_down_based_on_coverage(
					amr_name, coverage_annotation_list[i], summary_file,
					ref_up_info_lists[i], ref_down_info_lists[i], ref_amr_info_lists[i],
					params.assembler)
			average_precision+=precision
			average_sensitivity+=sensitivity
			logging.info('For "'+amr_name+'": sensitivity= '+str(sensitivity)+' precision = '+ str(precision))
		elif not coverage_annotation_list:
			ref_up_info_list, ref_amr_info_list, ref_down_info_list =\
				read_ref_annotations_from_db(amr_groups, amr_name)
			annotate_dir = params.output_dir+ANNOTATION_DIR+'/'+ANNOTATION_DIR+'_'+\
				str(params.seq_length)+'/annotation_'+restricted_amr_name+'_'+str(params.seq_length)
			if params.coverage_thr >=0:
				coverage_annotation = annotate_dir+'/coverage_annotation_'+\
					str(params.coverage_thr)+'_'+restricted_amr_name+'.csv'
			else:
				coverage_annotation = annotate_dir+'/trimmed_annotation_info_'+\
					restricted_amr_name+'.csv'
			sensitivity, precision = evaluate_sequences_up_down_based_on_coverage(
					amr_name, coverage_annotation, summary_file,
					ref_up_info_list, ref_down_info_list, ref_amr_info_list,
					params.assembler)
			average_precision+=precision
			average_sensitivity+=sensitivity
			logging.info('For "'+amr_name+'": sensitivity= '+str(sensitivity)+' precision = '+ str(precision))
	#Go over the members in overlap list
	overlap_amr_counter = 0
	for i, head in enumerate(heads):
		restricted_head = restricted_amr_name_from_modified_name(head)
		annotate_dir = params.output_dir+ANNOTATION_DIR+'/'+ANNOTATION_DIR+'_'+\
			str(params.seq_length)+'/annotation_'+restricted_head+'_'+str(params.seq_length)
		if params.coverage_thr >=0:
			coverage_annotation = annotate_dir+'/coverage_annotation_'+\
				str(params.coverage_thr)+'_'+restricted_head+'.csv'
		else:
			coverage_annotation = annotate_dir+'/trimmed_annotation_info_'+\
				restricted_head+'.csv'
		for amr_name in member_lists[i]:
			ref_up_info_list, ref_amr_info_list, ref_down_info_list =\
				read_ref_annotations_from_db(amr_groups, amr_name)
			sensitivity, precision = evaluate_sequences_up_down_based_on_coverage(
					amr_name, coverage_annotation, summary_file,
					ref_up_info_list, ref_down_info_list, ref_amr_info_list,
					params.assembler)
			average_precision+=precision
			average_sensitivity+=sensitivity
			logging.info('For "'+amr_name+'": sensitivity= '+str(sensitivity)+' precision = '+ str(precision))
			overlap_amr_counter+=1
	#go over the list of AMRs not found in the graph
	for amr_name in not_found_amr_names:
		ref_up_info_list, ref_amr_info_list, ref_down_info_list =\
			read_ref_annotations_from_db(amr_groups, amr_name)
		evaluate_sequences_up_down_based_on_coverage(amr_name, "", summary_file,
				ref_up_info_list, ref_down_info_list, ref_amr_info_list,
				params.assembler)

	return average_precision/(len(amr_files)+overlap_amr_counter+len(not_found_amr_names)),\
			average_sensitivity/(len(amr_files)+overlap_amr_counter+len(not_found_amr_names))

def construct_combined_annotation_file(coverage_annotation, amr_file, params):
	"""
	To create an annotation file containing both the annotations of sequences extracted
	from ref genomes as well as the ones  extracted from the graph --> this is used for
	visualizing the annotations
	Parameters:
		coverage_annotation: the file containing the annotation of sequences extracted from the graph
		amr_file: the file containing the AMR
		params: the list of parameters imported from params.py
	Return:
		a csv file containing both ref and graph annotations
	"""
	ref_seq_file = params.output_dir+AMR_DIR_NAME+'AMR_ref_neighborhood.fasta'
	_, amr_name = retrieve_AMR(amr_file)
	restricted_amr_name = restricted_amr_name_from_modified_name(amr_name)
	#initializ visual_annotation file
	visual_annotation_csv =os.path.dirname(coverage_annotation)+'/vis_annotation.csv'
	fd = open(visual_annotation_csv, 'w')
	vis_writer = csv.writer(fd)
	vis_writer.writerow(['seq_name', 'seq_value', 'seq_length', 'gene', 'coverage',
					'length', 'start_pos', 'end_pos', 'target_amr'])
	#Reading corresponding sequences in ref genomes
	ref_seq_list = []
	for record in SeqIO.parse(open(ref_seq_file,'r'),'fasta'):
		if record.id.startswith(amr_name+'::'):
			ref_seq_list.append(str(record.seq))
	#import pdb; pdb.set_trace()
	#annotate ref sequences
	annotate_dir = params.output_dir+'tmp_ref_annotations/'
	if not os.path.exists(annotate_dir):
		os.makedirs(annotate_dir)
	for index, ref_seq in enumerate(ref_seq_list):
		annotation_prefix = 'ref_'+ restricted_amr_name +'__'+str(index)
		seq_info = annotate_sequence(ref_seq+"\n", annotation_prefix, annotate_dir,
                    params.PROKKA_COMMAND_PREFIX, params.use_RGI, params.RGI_include_loose,
					delete_prokka_dir = True)
		#Add annotation of ref genomes
		for gene_info in seq_info:
			vis_writer.writerow(['ref_'+str(index+1), gene_info['seq_value'],
					len(gene_info['seq_value']), gene_info['gene'], gene_info['coverage'],
					gene_info['length'], gene_info['start_pos'],
					gene_info['end_pos'], gene_info['target_amr']])

	#add annotation of extracted sequences from graph
	with open(coverage_annotation) as fr:
		myreader = DictReader(fr)
		for row in myreader:
			vis_writer.writerow([row['seq_name'], row['seq_value'], row['seq_length'],
					row['gene'], row['coverage'], row['length'], row['start_pos'],
					row['end_pos'], row['target_amr']])

	fd.close()

	return visual_annotation_csv

def seq_annotation_trim_main(params, amr_files, all_seq_info_lists,
								annotation_files, visualize = False):
	"""
	The core function to filter and remove genes that their coverage difference
	from AMR coverage is above a threshold
	Prameters:
		params: the list of parameters imported from params.py
		amr_files: the list of files containing AMRs
		all_seq_info_lists: the list of annotations of neighborhood sequences extracted from the graph
		annotation_files: the files containing annotation info
		visualize: if True, visualize
	"""
	coverage_annotation_list = []
	for i, amr_file in enumerate(amr_files):
		restricted_amr_name = extract_name_from_file_name(amr_file)
		#remove some extracted sequences based on coverage consistency
		annotate_dir = params.output_dir+ANNOTATION_DIR+'/'+ANNOTATION_DIR+'_'+\
			str(params.seq_length)+'/annotation_'+restricted_amr_name+'_'+str(params.seq_length)
		coverage_annotation = ''
		remained_seqs = []
		if params.coverage_thr>0:
			coverage_annotation, remained_seqs = check_coverage_consistency_remove_rest_seq(\
								all_seq_info_lists[i],
								params.coverage_thr, restricted_amr_name, annotate_dir+'/')
		if visualize:
			# create an image presenting the annotations for all sequences
			if coverage_annotation!='':
				visual_annotation_csv = construct_combined_annotation_file(
						coverage_annotation, amr_file, params)
			else:
				visual_annotation_csv = construct_combined_annotation_file(
						annotation_files[i], amr_file, params)
			visual_annotation = annotate_dir+'/gene_comparison_'+str(params.coverage_thr)+'_'+restricted_amr_name+'.png'
			visualize_annotation(visual_annotation_csv, output=visual_annotation)
		if params.coverage_thr>0:
			coverage_annotation_list.append(coverage_annotation)
	return coverage_annotation_list

def seq_annotation_main(params, seq_files, path_info_files, amr_files):
	"""
	The core function for annotation of neighborhood sequences of all AMRs
	Parameters:
		params: the list of parameters extracted from params.py
		seq_files: the list of sequence files
		path_info_files: the list of files containing node info for all sequences
		amr_files: the list of files containing AMRs
	Return:

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
	#extract ref neighborhood annotation from the file
	if params.ref_genomes_available:
		df = pd.read_csv(params.ref_ng_annotations_file, skipinitialspace=True,  keep_default_na=False)
		amr_groups = df.groupby('target_amr')

	ref_up_info_lists = []
	ref_down_info_lists = []
	ref_amr_info_lists = []
	all_seq_info_lists = []
	annotation_files = []
	for amr_file in amr_files:
		restricted_amr_name = extract_name_from_file_name(amr_file)
		_, amr_name = retrieve_AMR(amr_file)
		ref_up_info_list = []
		ref_down_info_list = []
		ref_amr_info_list = []
		if params.ref_genomes_available:
			ref_up_info_list, ref_amr_info_list, ref_down_info_list =\
				read_ref_annotations_from_db(amr_groups, amr_name)
			ref_up_info_lists.append(ref_up_info_list)
			ref_down_info_lists.append(ref_down_info_list)
			ref_amr_info_lists.append(ref_amr_info_list)
		neighborhood_file, nodes_info_file = find_corrsponding_seq_path_file(restricted_amr_name,
								neighborhood_files, nodes_info_files, params.seq_length)
		if neighborhood_file == -1:
			logging.error('no sequence file for '+ amr_file +' was found! We looked for a file like '+restricted_amr_name)
			import pdb; pdb.set_trace()
			sys.exit()
		if params.multi_processor:
			all_seq_info_list, annotation_file=\
				neighborhood_annotation_parallel(amr_name, neighborhood_file,
					nodes_info_file, params.seq_length,
					ref_up_info_list, ref_down_info_list,
					params.output_dir, params.PROKKA_COMMAND_PREFIX,params.use_RGI,
					params.RGI_include_loose, '_'+restricted_amr_name,
					params.amr_identity_threshold, params.ref_genomes_available,
					params.core_num)
		else:
			all_seq_info_list, annotation_file =\
				neighborhood_annotation(amr_name, neighborhood_file,
					nodes_info_file, params.seq_length,
					ref_up_info_list, ref_down_info_list,
					params.output_dir, params.PROKKA_COMMAND_PREFIX,params.use_RGI,
					params.RGI_include_loose, '_'+restricted_amr_name,
					params.amr_identity_threshold, params.ref_genomes_available)
		all_seq_info_lists.append(all_seq_info_list)
		annotation_files.append(annotation_file)

	return all_seq_info_lists, ref_up_info_lists,\
			ref_down_info_lists, ref_amr_info_lists, annotation_files

def seq_evaluation_main(params, seq_files, amr_files, not_found_amr_names):
	"""
	The core function to compare sequences rather than their annotation
	Parameters:
		params: the list of parameters imported from the params.py
		seq_files: The list of sequence files
		amr_files: the list of files containing the AMRs
		not_found_amr_names: the list of AMRs not found in the graph which are available in the ref genomes
	Return:
		average precision and average sensitivity
	"""
	seq_identity_thr = 90
	logging.info("Sequence Evaluation ...")
	if seq_files:
		neighborhood_files = seq_files
	else:
		neighborhood_files = extract_files(params.ng_seq_files, 'please provide the \
			address of the files containing all extracted sequences from AMR neighborhood \
			in the assembly graph')
	#extract ref neighborhood sequences
	if not params.ref_genomes_available:
		logging.error('not able to evaluate as ref genome information is not available!')
		sys.exit()
	#Reading corresponding sequences in ref genomes
	ref_seq_file = params.output_dir+AMR_DIR_NAME+'AMR_ref_neighborhood.fasta'
	ref_up_seq_list = collections.defaultdict(list)
	ref_down_seq_list = collections.defaultdict(list)
	for record in SeqIO.parse(open(ref_seq_file,'r'),'fasta'):
		amr_name = str(record.id).split('::')[0]
		ref_up_seq, ref_down_seq = split_up_down_seq(str(record.seq))
		if ref_up_seq!='' and ref_up_seq not in ref_up_seq_list[amr_name]:
			ref_up_seq_list[amr_name].append(ref_up_seq)
		if ref_down_seq!='' and ref_down_seq not in ref_down_seq_list[amr_name]:
			ref_down_seq_list[amr_name].append(ref_down_seq)
	#initializ required file
	evaluation_csv =params.output_dir+'/sequence_evaluation.csv'
	with open(evaluation_csv, 'a') as ed:
		eval_writer = csv.writer(ed)
		eval_writer.writerow(['amr_name', 'type', 'extracted_seq', 'ref_seq', 'identity', 'coverage'])
	summary_file = params.output_dir+'summaryMetrics_sequence_'+str(seq_identity_thr)+'.csv'
	with open(summary_file,'a') as fd:
		writer = csv.writer(fd)
		writer.writerow(['AMR', 'Unique_TP#', 'FP#', 'Unique_True#', 'found#','sensitivity', 'precision'])
	#extract the list of AMR groups
	overlap_file_name = params.output_dir+AMR_DIR_NAME+AMR_OVERLAP_FILE
	heads, member_lists, _ = extract_info_from_overlap_file(overlap_file_name)
	#Going over amr files
	average_precision = 0
	average_sensitivity = 0
	for amr_file in amr_files:
		restricted_amr_name = extract_name_from_file_name(amr_file)
		_, amr_name = retrieve_AMR(amr_file)
		neighborhood_seq_file, _ = find_corrsponding_seq_path_file(restricted_amr_name,
								neighborhood_files, [], params.seq_length)
		#compare and evaluate
		precision, sensitivity = evaluate_pure_sequences(amr_name, neighborhood_seq_file,
					ref_up_seq_list[amr_name], ref_down_seq_list[amr_name],
					evaluation_csv, summary_file, seq_identity_thr)
		average_precision+=precision
		average_sensitivity+=sensitivity
		logging.info('For "'+amr_name+'": sensitivity= '+str(sensitivity)+' precision = '+ str(precision))
	#Going over the members in overlap list
	overlap_amr_counter = 0
	for i, head in enumerate(heads):
		restricted_head = restricted_amr_name_from_modified_name(head)
		neighborhood_seq_file, _ = find_corrsponding_seq_path_file(restricted_head,
								neighborhood_files, [], params.seq_length)
		for amr_name in member_lists[i]:
			#compare and evaluate
			precision, sensitivity = evaluate_pure_sequences(amr_name, neighborhood_seq_file,
						ref_up_seq_list[head], ref_down_seq_list[head],
						evaluation_csv, summary_file, seq_identity_thr)
			average_precision+=precision
			average_sensitivity+=sensitivity
			logging.info('For "'+amr_name+'": sensitivity= '+str(sensitivity)+' precision = '+ str(precision))
			overlap_amr_counter+=1
	#Going over the list of AMRs not found in the graph
	for amr_name in not_found_amr_names:
		precision, sensitivity = evaluate_pure_sequences(amr_name, '',
					ref_up_seq_list[amr_name], ref_down_seq_list[amr_name],
					evaluation_csv, summary_file, seq_identity_thr)

	average_precision = average_precision/(len(amr_files)+overlap_amr_counter+len(not_found_amr_names))
	average_sensitivity =average_sensitivity/(len(amr_files)+overlap_amr_counter+len(not_found_amr_names))
	logging.info('average_precision:'+str(average_precision)+' average_sensitivity:'+str(average_sensitivity))

	return average_precision, average_sensitivity

def sequence_neighborhood_main(params, gfa_file, graph_file, amr_seq_align_info):
	"""
	The core function to extract the neighborhood of AMRs
	Parameters:
		params: the list pf parameters imported from params.py
		gfa_file: the file containing the assembly graph
		graph_file: the file containing the assembly graph generated by this code
		amr_seq_align_info: the alignment info (AMR alignment in the graph)
	Return:
		the list of files containing extracted sequence and the details of nodes representing them
	"""
	seq_files = []
	path_info_files = []
	logging.info("Extracting neighborhood sequences with length = %s", params.seq_length)
	if not gfa_file:
		gfa_file = verify_file_existence(graph_file, params.gfa_file, \
				'please provide the address of the file containing the assembly graph')
	#remove paths from GFA file
	delete_lines_started_with('P', gfa_file)
	sequence_dir = params.output_dir+SEQ_DIR_NAME+'/'+SEQ_DIR_NAME+'_'+str(params.seq_length)+'/'
	if not os.path.exists(sequence_dir):
		os.makedirs(sequence_dir)
	if params.multi_processor:
		p_extraction = partial(neighborhood_sequence_extraction, gfa_file, params.seq_length,
							sequence_dir, params.BANDAGE_PATH,
							params.amr_identity_threshold, SEQ_NAME_PREFIX,
							params.path_node_threshold , params.path_seq_len_percent_threshold,
							params.max_kmer_size, params.assembler)
		with Pool(params.core_num) as p:
			lists = p.map(p_extraction, amr_seq_align_info)
		seq_files, path_info_files = zip(*lists)
	else:
		for amr_file in amr_seq_align_info:
			seq_file, path_info_file = neighborhood_sequence_extraction(gfa_file, params.seq_length,
								sequence_dir, params.BANDAGE_PATH,
								params.amr_identity_threshold, SEQ_NAME_PREFIX,
								params.path_node_threshold , params.path_seq_len_percent_threshold,
								params.max_kmer_size, params.assembler, amr_file)
			if seq_file:
				path_info_files.append(path_info_file)
				seq_files.append(seq_file)

	return seq_files, path_info_files

def extract_amr_info(params, graph_file, ref_amr_files):
	"""
	To read AMR info if available or generate them otherwise, including alignment of AMRs in the graph!
	Parameters:
		params: the list of parameters imported from params.py
		graph_file: the file containing the assembly graph
		ref_amr_file: the list of files containing AMRs found in the ref genomes
	Return:
		the list of unique AMR files (the heads of the groups) and their alignment info
		as well as the list of AMRs not found in the graph
	"""
	not_found_amr_names = []
	if params.ref_genomes_available:
		logging.info("Checking if ref AMR genes are available in the assembly graph ...")
		gfa_file = verify_file_existence(graph_file, params.gfa_file, \
			'please provide the address of the file containing the assembly graph')
		#if overlap.txt and alignment files have already generated just use them
		align_dir = params.output_dir+AMR_DIR_NAME+AMR_ALIGN_DIR
		overlap_file_name = params.output_dir+AMR_DIR_NAME+AMR_OVERLAP_FILE
		if os.path.exists(align_dir) and os.path.isfile(overlap_file_name):
			all_align_files = extract_files(align_dir, "the directory "+align_dir+" doesn't exist!")
			unique_amr_files, not_found_amr_names, amr_count = read_info_from_overlap_ref_files(
					overlap_file_name, ref_amr_files)
			#extract path_info for unique AMRs
			unique_amr_path_list = extract_path_info_for_amrs(all_align_files, unique_amr_files,
														amr_count, params.amr_identity_threshold)
		else:
			#create the alignment files for the AMRs
			if params.multi_processor:
				unique_amr_files, not_found_amr_names, unique_amr_path_list =\
				 				find_amrs_alignment_in_graph_parallel(gfa_file,
								params.output_dir, ref_amr_files, params.BANDAGE_PATH,
								params.amr_identity_threshold, params.core_num)
			else:
				unique_amr_files, not_found_amr_names, unique_amr_path_list =\
								find_amrs_alignment_in_graph(gfa_file,
								params.output_dir, ref_amr_files, params.BANDAGE_PATH,
								params.amr_identity_threshold)
	#if no ref is available detect AMRs in the assembly graph itself
	elif params.find_amr_genes:
		logging.info("Finding AMR genes in the assembly graph ...")
		gfa_file = verify_file_existence(graph_file, params.gfa_file, \
				'please provide the address of the file containing the assembly graph')
		if params.multi_processor:
			unique_amr_files, unique_amr_path_list =\
							find_all_amr_in_graph_parallel(gfa_file, params.output_dir,
							params.CARD_AMR_SEQUENCES, params.BANDAGE_PATH,
							params.amr_identity_threshold, params.core_num)
		else:
			unique_amr_files, unique_amr_path_list =\
							find_all_amr_in_graph(gfa_file, params.output_dir,
							params.CARD_AMR_SEQUENCES, params.BANDAGE_PATH,
							params.amr_identity_threshold)
	#if ref genomes are not available but the list of AMRs has already been found in the graph and is available
	else:
		align_dir = params.output_dir+AMR_DIR_NAME+AMR_ALIGN_DIR
		overlap_file_name = params.output_dir+AMR_DIR_NAME+AMR_OVERLAP_FILE
		#We assume that only unique AMRs (heads) are stored
		unique_amr_files = extract_files(params.amr_files, 'please provide the address of the AMR gene(s)')
		all_align_files = extract_files(align_dir, 'the alignments are not available!')
		heads, member_lists, unique_amr_list = extract_info_from_overlap_file(overlap_file_name)
		if len(unique_amr_list+heads)!=len(unique_amr_files):
			logging.error("inconsisteny between the list of provided unique AMRs and overlapfile info!")
			import pdb;pdb.set_trace()
		amr_count = len(heads) + len(unique_amr_list) + sum(len(list) for list in member_lists)
		#extract path_info for unique AMRs
		unique_amr_path_list = extract_path_info_for_amrs(all_align_files, unique_amr_files,
												amr_count, params.amr_identity_threshold)

	return unique_amr_files, not_found_amr_names, unique_amr_path_list

def full_pipeline_main(params):
	logging.info("Startting the pipeline ...")
	#Validate task values
	task_list = validate_task_values(params.task)
	# if params.artificial_amr_insertion and params.find_amr_genes:
	# 	logging.error("variables 'artificial_amr_insertion' and 'find_amr_genes' cannot be True at the same run!")
	# 	import pdb; pdb.set_trace()
	# 	sys.exit()
	graph_file =""
	metagenome_file = ""
	genome_amr_files = []
	path_info_files = []
	seq_files = []
	reads = ""
	gfa_file = None
	ref_amr_files = []
	ref_genome_files = []
	if params.ref_genomes_available and params.artificial_amr_insertion:
		if params.inserted_amr_file=='' or not os.path.isfile(params.inserted_amr_file):
			logging.error("ERROR: no AMR file has been provided!")
			import pdb; pdb.set_trace()
			sys.exit()

	if params.ref_genomes_available and (Pipeline_tasks.metagenome_creation.value in task_list or
										 not os.path.exists(params.output_dir+AMR_DIR_NAME)):
		ref_genome_files = extract_files(params.ref_genome_files, 'please provide the address of genome files')
		if Pipeline_tasks.metagenome_creation.value in task_list:
			logging.info("Creating the metagenome sample ...")
			if params.artificial_amr_insertion:
				genome_amr_files, metagenome_file = create_metagenome_with_amr_insertion(ref_genome_files,
							params.number_of_insertions, params.insertion_type, params.insertion_locations,
							params.inserted_amr_file, params.output_dir)
			else:
				metagenome_file = concatenate_files(ref_genome_files, params.metagenome_file)
		#if ref_genome is available and the directory containing detected AMR sequences
		#is not available then call find_amrs_in sample.py first to detect the AMRs from the ref sample
		if not os.path.exists(params.output_dir+AMR_DIR_NAME):
			#os.makedirs(params.output_dir+AMR_DIR_NAME)
			if find_annotate_amrs_in_ref(params.CARD_AMR_SEQUENCES, params.metagenome_file,
											params.PROKKA_COMMAND_PREFIX, params.use_RGI,
											params.RGI_include_loose, params.output_dir)==-1:
				print('please enter the path for the sample file and amr sequences file')
				import pdb; pdb.set_trace()
				sys.exit()
	if params.ref_genomes_available:
		ref_amr_files = extract_files(params.amr_files, 'please provide the address of the AMR gene(s)')

	if Pipeline_tasks.read_simulation.value in task_list:
		logging.info("Simulating reads ...")
		meta_file = verify_file_existence(metagenome_file, params.metagenome_file, \
					'please provide the address of metagenome file')
		reads = simulate_reads(meta_file, params.read_length, params.ART_PATH)

	if Pipeline_tasks.assembly.value in task_list:
		logging.info("Assembly ...")
		read_files = verify_file_existence(reads, params.reads, \
			'please provide the address of paired end reads')
		spade_output = params.output_dir + params.assembler_output_dir
		graph_file, contigs_file = do_assembly(read_files, params.SPADES_PATH,
										spade_output, params.spades_thread_num,
										params.spades_error_correction)
	#extract AMR and alignment information
	unique_amr_files, not_found_amr_names, unique_amr_path_list =\
			extract_amr_info(params, graph_file, ref_amr_files)
	if not unique_amr_files:
		logging.info("No AMR gene was found!")
		import pdb; pdb.set_trace()
		sys.exit()
	send_amr_align_info = False
	if unique_amr_path_list and len(unique_amr_path_list)==len(unique_amr_files):
		send_amr_align_info = True
	else:
		logging.error("AMR alignment info are not available")
		import pdb; pdb.set_trace()

	#create pairs of seq and align info
	amr_seq_align_info = []
	for i, amr_file in enumerate(unique_amr_files):
		amr_seq_align_info.append((amr_file, unique_amr_path_list[i]))
	# if Pipeline_tasks.graph_neighborhood.value in task_list:
	# 	graph_extraction_main(params, gfa_file, graph_file, amr_seq_align_files)

	if Pipeline_tasks.sequence_neighborhood.value in task_list:
		if MULTIPLE_SEQ_LENGTH:
			for seq_len in [500, 1000, 1500, 2000, 2500, 3000, 3500, 4000, 4500, 5000]:
				params.seq_length = seq_len
				seq_files, path_info_files = sequence_neighborhood_main(params, gfa_file,
						graph_file, amr_seq_align_info)
		else:
			seq_files, path_info_files = sequence_neighborhood_main(params, gfa_file,
					graph_file, amr_seq_align_info)

	if SEQ_EVAL:
		seq_evaluation_main(params, seq_files, unique_amr_files, not_found_amr_names)
		logging.info("All Done!")
		sys.exit()

	coverage_annotation_list = []
	ref_up_info_lists = []
	ref_down_info_lists = []
	ref_amr_info_lists = []
	if Pipeline_tasks.neighborhood_annotation.value in task_list:
		all_seq_info_lists, ref_up_info_lists,\
		ref_down_info_lists, ref_amr_info_lists, annotation_file_list =\
			seq_annotation_main(params, seq_files, path_info_files, unique_amr_files)
		#if MULTIPLE_COV_THR:
		if isinstance(params.coverage_thr, list) and len(coverage_thr)>1:
			coverage_thr_list = params.coverage_thr
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
				# for cov_thr in range(1, 50):
				# 	params.coverage_thr = cov_thr
				for cov_thr in coverage_thr_list:
					params.coverage_thr = cov_thr
					coverage_annotation_list = seq_annotation_trim_main(params, unique_amr_files,\
					 		all_seq_info_lists, annotation_file_list, False)
					if params.ref_genomes_available and Pipeline_tasks.neighborhood_evaluation.value in task_list:
						precision, sensitivity = annotation_evaluation_main(params, unique_amr_files,
									coverage_annotation_list, ref_up_info_lists,
									ref_down_info_lists, ref_amr_info_lists, not_found_amr_names)
						data.append([cov_thr, precision, 'Precision'])
						data.append([cov_thr, sensitivity, 'Sensitivity'])
						writer.writerow([cov_thr, precision, sensitivity])
			df = pd.DataFrame(data, columns = ['cov_thr', 'value', 'type'])
			sns.scatterplot(data=df, x = 'cov_thr', y = 'value', hue='type', style='type')
			plt.show()
		else:
			if isinstance(params.coverage_thr, list):
				params.coverage_thr = params.coverage_thr[0]
			coverage_annotation_list = seq_annotation_trim_main(params, unique_amr_files,\
				all_seq_info_lists, annotation_file_list, True)

	if params.ref_genomes_available and not isinstance(params.coverage_thr, list)\
		and Pipeline_tasks.neighborhood_evaluation.value in task_list:
		annotation_evaluation_main(params, unique_amr_files, coverage_annotation_list,
					ref_up_info_lists, ref_down_info_lists, ref_amr_info_lists,
					not_found_amr_names)

	logging.info("All Done!")

def create_arguments(params, parser):
	"""
	To create all required qrguments
	"""
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
	parser.add_argument('--main_dir', '-m', type = str, default=params.main_dir,
		help = 'the main dir to retrieve required files')
	parser.add_argument('--output_dir', '-O', type = str, default=params.output_dir,
		help = 'the output dir to store the results')
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
	parser.add_argument('--path_seq_len_percent_threshold', type = int, default = params.path_seq_len_percent_threshold,
		help = 'the threshold used for recursive pre_seq and post_seq until we have this percentage of the required length\
		 after which we just extract from the longest neighbor')
	parser.add_argument('--ref_genomes_available', type = str2bool, default = params.ref_genomes_available,
		help = 'Whether we have access to reference genome(s)')
	parser.add_argument('--multi_processor', type = str2bool, default = params.multi_processor,
		help = 'Whether to use multi processors for parallel programming')
	parser.add_argument('--core_num', type = int, default = params.core_num,
		help = 'the number of cores used in case of parallel programming')
	parser.add_argument('--coverage_thr', nargs="+", default = params.coverage_thr,
		help = 'coverage threshold to check if an annotated gene is truly AMR neighbor or just a false positive')
	parser.add_argument('--ref_ng_annotations_file', type = str, default = params.ref_ng_annotations_file,
		help = 'the file containing the annotation of all neighborhoods extracted from ref genomes.')
	parser.add_argument('--prokka_prefix', type = str, default = params.PROKKA_COMMAND_PREFIX,
		help = 'Set only if prokka is run through docker')

	return parser

def modify_params(params, args):
	"""
	"""
	#updating params values
	params.task = args.task
	params.amr_files = args.amr_files
	params.ref_genome_files = args.ref_genome_files
	params.main_dir = args.main_dir
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
	params.reads = args.reads
	params.spades_error_correction = args.spades_error_correction
	params.use_RGI = args.use_RGI
	params.RGI_include_loose = args.RGI_include_loose
	params.find_amr_genes = args.find_amr_genes
	params.metagenome_file =  args.metagenome_file
	params.artificial_amr_insertion = args.artificial_amr_insertion
	params.amr_identity_threshold = args.amr_identity_threshold
	params.path_node_threshold = args.path_node_threshold
	params.path_seq_len_percent_threshold = args.path_seq_len_percent_threshold
	params.ref_genomes_available = args.ref_genomes_available
	params.coverage_thr = args.coverage_thr
	params.multi_processor = args.multi_processor
	params.core_num = args.core_num
	params.ref_ng_annotations_file = args.ref_ng_annotations_file
	params.PROKKA_COMMAND_PREFIX = args.prokka_prefix

	return params

if __name__=="__main__":
	import params
	text = 'This code is used to find the context of a given AMR gene'
	parser = argparse.ArgumentParser(description=text)
	parser = create_arguments(params, parser)
	args = parser.parse_args()
	params = modify_params(params, args)
	log_name = 'logger_'+datetime.datetime.now().strftime('%Y-%m-%d_%H-%M')+'.log'
	initialize_logger(params.main_dir, log_name)
	logging.info(str(params.__dict__))
	#create the output directory; if it exists, delete it and create a new one
	if not os.path.exists(params.output_dir):
		os.makedirs(params.output_dir)
	else:
		try:
			shutil.rmtree(params.output_dir)
		except OSError as e:
			logging.error("Error: %s - %s." % (e.filename, e.strerror))
		os.makedirs(params.output_dir)
	full_pipeline_main(params)
