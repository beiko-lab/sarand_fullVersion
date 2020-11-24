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
import gfapy
import re
import argparse
import difflib
import datetime
import csv
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

from extract_neighborhood import neighborhood_graph_extraction, neighborhood_sequence_extraction,\
							extract_amr_neighborhood_in_ref_genome, extract_nodes_in_path
# from extract_neighborhood_3_beforeBothDirRecursively import neighborhood_graph_extraction,\
#  							neighborhood_sequence_extraction,\
# 							extract_amr_neighborhood_in_ref_genome, extract_nodes_in_path
from find_seq_in_contigs import find_sequence_match
from annotation_visualization import visualize_annotation

"""
# The list of valid tasks
# to use them set task in params.py:
# for a single task: insert its value as the single item in task list
# for a range of tasks: insert the value of start and end as the first and the second (last) item in the list
# for all tasks: task = [0]
"""
class Pipeline_tasks(enum.Enum):
	all = 0
	metagenome_creation = 1
	read_simulation = 2
	assembly = 3
	graph_neighborhood = 4
	sequence_neighborhood = 5
	neighborhood_comparison = 6

# Whether to insert the AMR sequence in arandom or pre-defined location of the genome
# in case of amr artificial insertion
class Insertion_type(enum.Enum):
	random = 1
	assigned = 2

ASSEMBLY_FILE = 'assembly_graph_with_scaffolds.gfa'
CONTIG_FILE = 'contigs.fasta'
ALL_AMR_SEQUENCES ='nucleotide_fasta_protein_homolog_model.fasta'
AMR_DIR_NAME = 'AMR/'
SUBGRAPH_DIR_NAME = 'subgraphs/'
SEQ_DIR_NAME = 'sequences_info/'
SEQ_NAME_PREFIX = 'ng_sequences_'

def initialize_logger(output_dir, file_name = 'logfile.log'):
	"""
	"""
	#logging_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)),'logs')
	logging_dir = os.path.join(output_dir,'logs')
	try:
		os.makedirs(logging_dir)
	except OSError:
		pass
	log_file_path = os.path.join(logging_dir, file_name)
	logging.basicConfig(level=logging.DEBUG)
	log_formatter = logging.Formatter(
		"%(asctime)s [%(levelname)-5.5s]	%(message)s"
	)
	root_logger = logging.getLogger()
	file_handler = logging.FileHandler(log_file_path)
	file_handler.setFormatter(log_formatter)
	root_logger.addHandler(file_handler)
	console_handler = logging.StreamHandler(sys.stdout)
	console_handler.setFormatter(log_formatter)
	root_logger.addHandler(console_handler)

def str2bool(v):
	"""
	To convert a string to a boolean value
	Parameter:
		v: the input string
	Return:
		the converted boolean value
	"""
	if isinstance(v, bool):
		return v
	if v.lower() in ('yes', 'true', 't', 'y', '1'):
		return True
	elif v.lower() in ('no', 'false', 'f', 'n', '0'):
		return False
	else:
		raise argparse.ArgumentTypeError('Boolean value expected.')

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
	"\nneighborhood_comparison = "+str(Pipeline_tasks.neighborhood_comparison.value)
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

def print_params(params):
	"""
	"""
	mylog = '\n'.join(param for param in params)
	logging.info(mylog)

def verify_file_existence(myfile, param_file, message):
	"""
	This code is used to check if the required file exits either through running
	previous steps of the pipeline or passed by user (or in the param list)
	Parameters:
		myfile:	the file was supposed to be generated via previous steps of the pipeline
		param_file:	the file set by the user
		message:	error message in case no instance of the required file exists
	Return:
		the valid instance of the parameter to work with
	"""
	if myfile!="":
		return myfile
	elif param_file !="":
		if os.path.isfile(param_file):
			return param_file
	logging.error("ERROR: "+message)
	#print("ERROR: "+message)
	sys.exit()



def retrieve_AMR(file_path):
	"""
	To read the AMR gene from the text file.
	Parameters:
		file_path:	the address of the file containing the AMR gene
	Return:
		the sequence of the AMR gene in lower case
	"""
	with open(file_path) as fp:
		for i, line in enumerate(fp):
			#skip comment line
			if line.startswith('>'):
				continue
			return line

def write_info_in_annotation_file(annotation_writer, visual_annotation_writer,
					seq_description, seq, contig_name, seq_info, use_RGI, found):
	"""
	To write annotation details into files
	Parameters:
		annotation_writer:	annotation file containing all annotations
		visual_annotation_writer: annotation file containing unique annotations
		seq_description: a small description of the sequence used for naming
		seq: the extracted dna sequence that has been annotated
		seq_info: annotation info
		contig_name: the name of contig matching this extracted sequence (if there is any contig)
		use_RGI: True if RGI has been used to annotate AMRs
		found: True if the annotation info has already found in other annotated sequences
	"""
	if use_RGI:
		annotation_writer.writerow([seq_description, seq, len(seq),
							contig_name, seq_info['gene'], seq_info['prokka_gene_name'],
							seq_info['product'], seq_info['length'],
							seq_info['start_pos'], seq_info['end_pos'],
							seq_info['RGI_prediction_type']])
		if not found:
			visual_annotation_writer.writerow([seq_description, seq, len(seq),
								contig_name, seq_info['gene'], seq_info['prokka_gene_name'],
								seq_info['product'], seq_info['length'],
								seq_info['start_pos'], seq_info['end_pos'],
								seq_info['RGI_prediction_type']])
	else:
		annotation_writer.writerow([seq_description, seq, len(seq),
							contig_name, seq_info['gene'],
							seq_info['product'], seq_info['length'],
							seq_info['start_pos'], seq_info['end_pos']])
		if not found:
			visual_annotation_writer.writerow([seq_description, seq, len(seq),
								contig_name, seq_info['gene'],
								seq_info['product'], seq_info['length'],
								seq_info['start_pos'], seq_info['end_pos']])


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
	return read1, read2

def do_assembly(read1, read2, spades_path, output_dir, thread_num = 16, error_correction = True):
	"""
	To call MetaSPAdes and run assembly
	Parameters:
		read1, read2:	the address of paired end reads
		spades_path:	the address of SPAdes to run from
		output_dir:		the path for the output directory
		thread_num:		the number of threads used in MetaSPAdes
		error_correction:to choose if Spades's error correction feature is used or not
	Return:
		the address of the file containing assembly graph (gfa_file) and
			contigs (contig_file)
	"""
	if error_correction:
		command = 'python ' + spades_path + ' -1 '+read1+' -2 '+read2+' --meta --threads ' \
					+ str(thread_num) + ' -o '+ output_dir
	else:
		command = 'python ' + spades_path + ' -1 '+read1+' -2 '+read2+' --meta \
		--only-assembler --threads ' + str(thread_num) + ' -o '+ output_dir
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
		output_file_name +" --input_type protein --clean"
	if include_loose:
		command+=" --include_loose --exclude_nudge"
	os.system(command)

	seq_info_list = []
	if os.path.isfile(output_file_name + '.txt'):
		with open(output_file_name + '.txt', newline = '') as rgi_file:
			rgi_reader = csv.reader(rgi_file, delimiter='\t')
			next(rgi_reader)
			for row in rgi_reader:
				seq_info = {'ORF_ID':row[0], 'gene':row[8].strip(),
				'prediction_type':row[5].strip(), 'best_identities':float(row[9])}
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
	seq_file_name = create_fasta_file(seq, output_dir)

	prokka_dir = output_dir+'prokka_dir_'+seq_description+'_'+datetime.datetime.now().strftime('%Y-%m-%d_%H-%M')
	if os.path.exists(prokka_dir):
		try:
			shutil.rmtree(prokka_dir)
		except OSError as e:
			logging.error("Error: %s - %s." % (e.filename, e.strerror))
			#print("Error: %s - %s." % (e.filename, e.strerror))

	command = prokka_prefix + 'prokka --metagenome --outdir '+\
		prokka_dir+' --prefix mygenome '+seq_file_name
	os.system(command)

	RGI_output_list = None
	if use_RGI:
		RGI_output_list = run_RGI(prokka_dir+"/mygenome.faa", output_dir, seq_description,
								RGI_include_loose)

	#Go over Prokka's output files and extract required information
	seq_info_list = []
	with open(prokka_dir+"/mygenome.tsv") as tsvfile:
		reader = csv.reader(tsvfile, delimiter='\t')
		#skip the header
		next(reader)
		for row in reader:
			seq_info = {'locus_tag':row[0].strip(), 'gene':row[3].strip(),'length':row[2].strip(),
						'product':row[6].strip(),'start_pos':None, 'end_pos':None,
						'prokka_gene_name':row[3].strip(), 'RGI_prediction_type':None}
			seq_info_list.append(seq_info)
	counter = 0
	with open(prokka_dir+"/mygenome.tbl", 'r') as read_obj:
		for line in read_obj:
			if line[0].isdigit():
				cells = line.split('\t')
				seq_info_list[counter]['start_pos'] = int(cells[0])
				seq_info_list[counter]['end_pos'] = int(cells[1])
				counter+=1

	#incorporate RGI findings into Prokka's
	if RGI_output_list:
		for item in RGI_output_list:
			for seq_info in seq_info_list:
				if item['ORF_ID'].split(' ')[0]==seq_info['locus_tag']:
					seq_info['gene'] = item['gene']
					seq_info['RGI_prediction_type'] = item['prediction_type']
					break

	# #remove temporary files and folder
	# os.remove(seq_file_name)
	# try:
	# 	shutil.rmtree(prokka_dir)
	# except OSError as e:
	# 	print("Error: %s - %s." % (e.filename, e.strerror))

	return seq_info_list

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
		if len(seq_list)==len(seq_info_list):
			identical_rows = 0
			for i, seq_info in enumerate(seq_info_list):
				if (seq_info['gene']==seq_list[i]['gene'] or\
					(seq_info['gene']=='' and seq_list[i]['gene']=='UNKNOWN')) and\
					seq_info['length']==seq_list[i]['length'] and\
				 	seq_info['start_pos']==seq_list[i]['start_pos'] and\
					seq_info['end_pos']==seq_list[i]['end_pos']:
					identical_rows+=1
			if identical_rows == len(seq_list):
				found = True
				break

	return found

def neighborhood_comparison(amr_file, ref_genome_files, neighborhood_seq_file, \
								contig_file, seq_length, output_dir, prokka_prefix,\
								use_RGI = True, RGI_include_loose = False, output_name ='',
								amr_threshold = 95):
	"""
	To annotate reference genomes (a piece extracted around the AMR gene) as well as
		extracted neighborhood sequences from assembly graph, summarize the results
		in a couple of formats, visualize them (and compare them : ???????to be implemented?????)
	Parameters:
		amr_file:	the address of the file containing the AMR sequence
		ref_genome_files:	the list of address of reference genomes after
			inserting the AMR gene in them
		neighborhood_seq_file:	the address of the file containing all extracted
		 	neighborhood sequences from assembly graph
		contig_file:	the address of file containing the contigs of assembly
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
			visual_annotation_info, gene_file_name, product_file_name, visual_annotation)
	"""
	# initializing required files and directories
	annotate_dir = output_dir+"annotation/annotation"+output_name+'_'+str(seq_length)
	if os.path.exists(annotate_dir):
		try:
			shutil.rmtree(annotate_dir)
		except OSError as e:
			logging.error("Error: %s - %s." % (e.filename, e.strerror))
			#print("Error: %s - %s." % (e.filename, e.strerror))
	os.makedirs(annotate_dir)
	annotation_detail_name = annotate_dir+'/annotation_detail'+output_name+'.csv'
	visual_annotation_info_name = annotate_dir+'/visual_annotation_info'+output_name+'.csv'
	annotation_detail = open(annotation_detail_name, mode='w', newline='')
	visual_annotation_info = open(visual_annotation_info_name, mode='w', newline='')
	#annotation_writer = csv.writer(annotation_detail, delimiter=',', quoting=csv.QUOTE_MINIMAL)
	annotation_writer = csv.writer(annotation_detail)
	visual_annotation_writer = csv.writer(visual_annotation_info)
	if use_RGI:
		annotation_writer.writerow(["seq_name", "seq_value", "seq_length",\
								"matched_contig", "gene", "prokka_gene_name", "product", \
								"length", "start_pos", "end_pos", 'RGI_prediction_type'])
		visual_annotation_writer.writerow(["seq_name", "seq_value", "seq_length",\
								"matched_contig", "gene", "prokka_gene_name", "product", \
								"length", "start_pos", "end_pos", 'RGI_prediction_type'])
	else:
		annotation_writer.writerow(["seq_name", "seq_value", "seq_length",\
								"matched_contig", "gene", "product", "length", "start_pos", "end_pos"])
		visual_annotation_writer.writerow(["seq_name", "seq_value", "seq_length",\
								"matched_contig", "gene", "product", "length", "start_pos", "end_pos"])
	gene_file_name = annotate_dir+'/seq_comparison_genes'+output_name+'.txt'
	gene_file = open(gene_file_name, 'w')
	product_file_name = annotate_dir+'/seq_comparison_products'+output_name+'.txt'
	product_file = open(product_file_name, 'w')

	#find neighborhood sequence in ref genome(s) and annotate using prokka and RGI
	amr_seq = retrieve_AMR(amr_file)
	ref_file_name = annotate_dir+'/ref_neighborhood_sequences'+output_name+'_'+\
		str(seq_length)+'_'+datetime.datetime.now().strftime('%Y-%m-%d_%H-%M')+'.txt'
	ref_file = open(ref_file_name, 'w')
	for genome_file in ref_genome_files:
		seq_list = extract_amr_neighborhood_in_ref_genome(amr_seq, genome_file, seq_length, amr_threshold)
		#AMR sequence might be in multiple places in ref genome so we have a list of
		#extracted neighborhood sequences instead of one sequence
		for index, seq in enumerate(seq_list):
			genome_name = os.path.splitext(os.path.basename(genome_file))[0]
			ref_file.write("> " + genome_name+'_'+str(index+1) + " :\n")
			ref_file.write(seq+"\n")
			seq_description = 'ref_'+genome_name+'_'+str(index+1)
			seq_info_list = annotate_sequence(seq+"\n", seq_description, annotate_dir+'/',
										prokka_prefix, use_RGI, RGI_include_loose)
			myLine1 = myLine2 = seq_description + ':\t'
			for seq_info in seq_info_list:
				write_info_in_annotation_file(annotation_writer, visual_annotation_writer,
											seq_description, seq, "", seq_info, use_RGI, False)
				if seq_info['gene']=='':
					seq_info['gene']='UNKNOWN'
				myLine1+=seq_info['gene']+'---'
				myLine2+=seq_info['product']+'---'
			gene_file.write(myLine1[:-3]+'\n')
			product_file.write(myLine2[:-3]+'\n')
	ref_file.close()
	logging.info("NOTE: The list of neighborhood sequences in reference genome(s) has\
	 		been stroed in " + ref_file_name)
	#print("NOTE: The list of neighborhood sequences in reference genome(s) has\
	 #		been stroed in " + ref_file_name)

	#annotate the sequences extraced from assembly graph
	counter = 1
	all_seq_info_lists =[]
	logging.info('Reading '+ neighborhood_seq_file + ' for '+ amr_file)
	#print('Reading '+ neighborhood_seq_file + ' for '+ amr_file )
	with open(neighborhood_seq_file, 'r') as read_obj:
		for line in read_obj:
			if line.startswith('>') or line.startswith('Path') or line.startswith('The'):
				continue
			#find out if there is a match for it in contig file
			contig_list = find_sequence_match(line, contig_file, annotate_dir, '_'+str(counter))
			contig_name = ""
			for contig in contig_list:
				if int(float(contig.identity))>=amr_threshold and int(float(contig.matched_length)/(len(line)-1)*100)>=amr_threshold:
					contig_name = contig.name
					break
			seq_description = 'extracted'+str(counter)
			seq_info_list = annotate_sequence(line, seq_description, annotate_dir+'/',
											prokka_prefix, use_RGI, RGI_include_loose)
			found = seq_annotation_already_exist(seq_info_list, all_seq_info_lists)
			if not found:
				all_seq_info_lists.append(seq_info_list)
			myLine1 = myLine2 = seq_description +':\t'
			for seq_info in seq_info_list:
				write_info_in_annotation_file(annotation_writer, visual_annotation_writer,
											seq_description, line[:-1], contig_name,
											seq_info, use_RGI, found)
				if seq_info['gene']=='':
					seq_info['gene']='UNKNOWN'
				myLine1+=seq_info['gene']+'---'
				myLine2+=seq_info['product']+'---'
			gene_file.write(myLine1[:-3]+'\n')
			product_file.write(myLine2[:-3]+'\n')
			counter+=1

	annotation_detail.close()
	visual_annotation_info.close()
	gene_file.close()
	product_file.close()
	logging.info("NOTE: The comparison of neighborhood sequences are available in " +\
	 		annotation_detail_name+", "+gene_file_name+", "+product_file_name)
	#print("NOTE: The comparison of neighborhood sequences are available in " +\
	# 		annotation_detail_name+", "+gene_file_name+", "+product_file_name)

	# create an image presenting the annotations for all sequences
	visual_annotation = annotate_dir+'/gene_comparison_'+output_name+'.png'
	if len(all_seq_info_lists)<20:
		visualize_annotation(visual_annotation_info_name, output=visual_annotation)

	return annotation_detail_name, gene_file_name, product_file_name, visual_annotation

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
	return found, paths_info

def amr_path_overlap(found_amr_paths, new_paths, overlap_percent = 95):
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
					percent = (1 - (float(diff_length)/(new_path['end_pos'] - new_path['start_pos'])))*100
					if percent >= overlap_percent:
						found = True
						id_list.append(i)
						break
			if found:
				break
	if len(id_list)==len(new_paths):
		return True, id_list
	return False, None

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
	align_dir = output_dir+AMR_DIR_NAME+'alignments/'
	if not os.path.exists(align_dir):
		os.makedirs(align_dir)

	amr_name = ''
	found_amr_seqs = []
	found_amr_names = []
	found_amr_paths = []
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
			found, paths_info = is_there_amr_in_graph(amr_file, amr_name_processed, gfa_file,
											align_dir, bandage_path, amr_threshold)
			if found:
				overlap, amr_ids =  amr_path_overlap(found_amr_paths, paths_info)
				if not overlap:
					found_amr_seqs.append(line)
					amr_info = {'name':amr_name, 'overlap_list':[]}
					found_amr_names.append(amr_info)
					found_amr_paths.append(paths_info)
				else:
					# add this AMR to the right group of AMRs all having overlaps
					for id in amr_ids:
						found_amr_names[id]['overlap_list'].append(amr_name_processed1)

	# write information (the sequence of found AMRs that don't have overlaped paths with others
	# + the list of groups in each all AMRs have overlaped paths) into files
	AMR_dir = output_dir+AMR_DIR_NAME+'sequences/'
	if not os.path.exists(AMR_dir):
		os.makedirs(AMR_dir)
	overlap_file_name = output_dir+AMR_DIR_NAME+'overlaps.txt'
	overlap_file = open(overlap_file_name, 'w')
	amr_files = []
	for i, seq in enumerate(found_amr_seqs):
		amr_name1 = found_amr_names[i]['name'].split('[')[0].split('|')[-1].strip()
		amr_name = ''.join(e for e in amr_name1 if e.isalpha() or e.isnumeric() or e=='_' or e=='-')
		amr_file = create_fasta_file(seq, AMR_dir, found_amr_names[i]['name'], amr_name)
		amr_files.append(amr_file)
		if found_amr_names[i]['overlap_list']:
			overlap_file.write(amr_name1+":\n")
			overlap_file.write(', '.join(e for e in found_amr_names[i]['overlap_list']))
			overlap_file.write("\n")
	overlap_file.close()

	return amr_files

def extract_files(gfiles, message):
	"""
	To extract file(s) address from an object
	# if gfiles is a list and the first item is a file address (it would be more
	# accurate to check this for all items) return gfiles
	# else if gfiles is a file address return [gfiles] as a list with one item
	# else if gfiles is a directory address return the list of files in it
	Parameters:
		gfiles:		a string or list of strings
		message:	an error message in case that no file was extracted successfully
	Return:
		the list of file(s) address
	"""
	if isinstance(gfiles, list):
		#check if these are files (not directories)
		if os.path.isfile(gfiles[0]):
			return gfiles
		else:
			logging.error(message)
			#print(message)
			sys.exit()
		#elif os.path.isdir(gfiles[0])
	elif os.path.isfile(gfiles):
		return [gfiles]
	elif os.path.isdir(gfiles):
		myfiles = [os.path.join(gfiles, f) for f in os.listdir(gfiles) \
							if os.path.isfile(os.path.join(gfiles, f))]
		return myfiles
	else:
		logging.error(message)
		#print(message)
		sys.exit()

def find_corrsponding_seq_file(amr_name, sequences_file_names, seq_length):
	"""
	To return the name of a sequence file (from sequences_file_names)
	dedicated to sequences extracted with a given length (seq_length)
	from a given amr sequence (amr_name)
	"""
	for file_name in sequences_file_names:
		if SEQ_NAME_PREFIX+amr_name+'_'+str(seq_length) in file_name:
			return file_name
	return -1

def main(params):
	logging.info("Startting the pipeline ...")
	#Validate task values
	task_list = validate_task_values(params.task)
	if params.artificial_amr_insertion and params.find_amr_genes:
		logging.error("variables 'artificial_amr_insertion' and 'find_amr_genes' cannot be True at the same run!")
		#print("variables 'artificial_amr_insertion' and 'find_amr_genes' cannot be True at the same run!")
		sys.exit()

	graph_file =""
	contigs_file = ""
	read1 = ""
	read2 = ""
	metagenome_file = ""
	genome_amr_files = []
	seq_files = []
	gfa_file = None
	amr_files = []
	ref_genome_files = extract_files(params.ref_genome_files, 'please provide the address of genome files')

	if not params.find_amr_genes:
		amr_files = extract_files(params.amr_files, 'please provide the address of the AMR gene(s)')
	if params.artificial_amr_insertion:
		if not amr_files:
			logging.error("ERROR: no AMR file has been provided!")
			#print("ERROR: no AMR file has been provided!")
			sys.exit()
		# for now, we only accept one AMR gene to be inserted in the genomes
		if len(amr_files)>1:
			logging.error("ERROR: for the artificial AMR insertion, we can't accept more than one AMR gene sequence!")
			#print("ERROR: for the artificial AMR insertion, we can't accept more than one AMR gene sequence!")
			sys.exit()

	if Pipeline_tasks.metagenome_creation.value in task_list:
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
		read1, read2 = simulate_reads(meta_file, params.read_length, params.ART_PATH)

	if Pipeline_tasks.assembly.value in task_list:
		logging.info("Assembly ...")
		read1_file = verify_file_existence(read1, params.read1, \
			'please provide the address of both paired end reads')
		read2_file = verify_file_existence(read2, params.read2, \
			'please provide the address of both paired end reads')
		spade_output = params.output_dir + params.spades_output_dir
		graph_file, contigs_file = do_assembly(read1_file, read2_file, params.SPADES_PATH,
										spade_output, params.spades_thread_num,
										params.spades_error_correction)

	if params.find_amr_genes:
		logging.info("Finding AMR genes in the assembly graph ...")
		gfa_file = verify_file_existence(graph_file, params.gfa_file, \
				'please provide the address of the file containing the assembly graph')
		amr_files = find_all_amr_in_graph(gfa_file, params.output_dir,
							params.PATH_PREFIX+ALL_AMR_SEQUENCES, params.BANDAGE_PATH,
							params.amr_identity_threshold)

	if Pipeline_tasks.graph_neighborhood.value in task_list:
		logging.info("Extracting neighborhood subgraphs ...")
		if not gfa_file:
			gfa_file = verify_file_existence(graph_file, params.gfa_file, \
					'please provide the address of the file containing the assembly graph')
		if not os.path.exists(params.output_dir+SUBGRAPH_DIR_NAME):
			os.makedirs(params.output_dir+SUBGRAPH_DIR_NAME)
		for amr_file in amr_files:
			output_name = 'ng_subgraph_'+os.path.splitext(os.path.basename(amr_file))[0]
			_ = neighborhood_graph_extraction(amr_file, gfa_file, params.graph_distance,
				 				params.output_dir+SUBGRAPH_DIR_NAME, params.BANDAGE_PATH,
								params.amr_identity_threshold, output_name)

	if Pipeline_tasks.sequence_neighborhood.value in task_list:
		logging.info("Extracting neighborhood sequences with length = %s", params.seq_length)
		if not gfa_file:
			gfa_file = verify_file_existence(graph_file, params.gfa_file, \
					'please provide the address of the file containing the assembly graph')
		if not os.path.exists(params.output_dir+SEQ_DIR_NAME):
			os.makedirs(params.output_dir+SEQ_DIR_NAME)
		for amr_file in amr_files:
			logging.info('amr_file = '+amr_file)
			#print('amr_file = '+amr_file)
			output_name = SEQ_NAME_PREFIX+os.path.splitext(os.path.basename(amr_file))[0]
			seq_file = neighborhood_sequence_extraction(amr_file, gfa_file, params.seq_length,
								params.output_dir+SEQ_DIR_NAME, params.BANDAGE_PATH,
								params.amr_identity_threshold, output_name,
								params.path_node_threshold , params.path_seq_len_percent_threshod)
			if seq_file:
				seq_files.append(seq_file)

	if Pipeline_tasks.neighborhood_comparison.value in task_list:
		logging.info("Neighborhood Comparison ...")
		if seq_files:
			neighborhood_files = seq_files
		else:
			neighborhood_files = extract_files(params.neighborhood_seq_files, 'please provide the \
				address of the files containing all extracted sequences from AMR neighborhood \
				in the assembly graph')
		contig_file = verify_file_existence(contigs_file, params.contig_file,
				'please provide the address of the file containing contigs after assembly')
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
		for amr_file in amr_files:
			amr_name = os.path.splitext(os.path.basename(amr_file))[0]
			neighborhood_file = find_corrsponding_seq_file(amr_name, neighborhood_files, params.seq_length)
			if neighborhood_file == -1:
				logging.error('no sequence file for the corresponding amr file was found!')
				#print('no sequence file for the corresponding amr file was found!')
				sys.exit()
			annotation_detail, gene_annotation, product_annotation, visual_annotation=\
				neighborhood_comparison(amr_file, genome_files,neighborhood_file,
					contig_file, params.seq_length, params.output_dir,params.PROKKA_COMMAND_PREFIX,
					params.use_RGI, params.RGI_include_loose, '_'+amr_name, params.amr_identity_threshold)

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
		"\nneighborhood_comparison = "+str(Pipeline_tasks.neighborhood_comparison.value))
	parser.add_argument('--amr_files','-A', type=str, default = params.amr_files,
		help = 'the path of the file(s) containing the AMR gene sequence(s)')
	parser.add_argument('--ref_genome_files', nargs="+", default=params.ref_genome_files,
		help = 'the ddress of reference genomes that AMR genome will be inserted in them')
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
	parser.add_argument('--spades_output_dir',type=str, default=params.spades_output_dir,
		help = 'the output dir to store MetaSPAdes results')
	parser.add_argument('--graph_distance', '-D', type = int, default=params.graph_distance,
		help = 'the maximum distance of neighborhood nodes to be extracted from the AMR gene')
	parser.add_argument('--seq_length', '-L', type = int, default=params.seq_length,
		help = 'the length of AMR gene\'s neighbourhood to be extracted')
	parser.add_argument('--neighborhood_seq_files', nargs="+", default = params.neighborhood_seq_files,
		help = 'the address of the files containing all extracted neighborhood sequences in assembly graph')
	parser.add_argument('--gfa_file', type = str, default = params.gfa_file,
		help = 'the address of the file for assembly graph')
	parser.add_argument('--contig_file', type = str, default = params.contig_file,
		help = 'the address of the file containing contigs after assembly')
	parser.add_argument('--genome_amr_files', nargs="+", default = params.genome_amr_files,
		help = 'the address of the files containing genome after AMR insertion')
	parser.add_argument('--read1', type = str, default = params.read1,
		help = 'the address of the files containing the first file in paired-end reads')
	parser.add_argument('--read2', type = str, default = params.read2,
		help = 'the address of the files containing the second file in paired-end reads')
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

	args = parser.parse_args()

	#updating params values
	#parameters = {}
	params.task = args.task
	params.amr_files = args.amr_files
	params.ref_genome_files = args.ref_genome_files
	params.output_dir = args.output_dir
	params.number_of_insertions = args.number_of_insertions
	params.insertion_type = args.insertion_type
	params.insertion_locations = args.insertion_locations
	params.read_length = args.read_length
	params.spades_thread_num = args.spades_thread_num
	params.spades_output_dir = args.spades_output_dir
	params.graph_distance = args.graph_distance
	params.seq_length = args.seq_length
	params.neighborhood_seq_files = args.neighborhood_seq_files
	params.gfa_file = args.gfa_file
	params.contig_file = args.contig_file
	params.genome_amr_files = args.genome_amr_files
	params.read1 = args.read1
	params.read2 = args.read2
	params.spades_error_correction = args.spades_error_correction
	params.use_RGI = args.use_RGI
	params.RGI_include_loose = args.RGI_include_loose
	params.find_amr_genes = args.find_amr_genes
	params.metagenome_file =  args.metagenome_file
	params.artificial_amr_insertion = args.artificial_amr_insertion
	params.amr_identity_threshold = args.amr_identity_threshold
	params.path_node_threshold = args.path_node_threshold
	params.path_seq_len_percent_threshod = args.path_seq_len_percent_threshod

	log_name = 'logger_'+datetime.datetime.now().strftime('%Y-%m-%d_%H-%M')+'.log'
	initialize_logger(params.output_dir, log_name)
	logging.info(str(params.__dict__))
	#print_params(params)

	main(params)
