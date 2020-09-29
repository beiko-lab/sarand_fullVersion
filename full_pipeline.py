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

import params
from extract_neighborhood import neighborhood_graph_extraction, neighborhood_sequence_extraction,\
							extract_amr_neighborhood_in_ref_genome
from find_seq_in_contigs import find_sequence_match
from annotation_visualization import visualize_annotation

class Pipeline_tasks(enum.Enum):
	all = 0
	metagenome_creation = 1
	read_simulation = 2
	assembly = 3
	graph_neighborhood = 4
	sequence_neighborhood = 5
	neighborhood_comparison = 6
	#contig_matching = 7

class Insertion_type(enum.Enum):
	random = 1
	assigned = 2

ASSEMBLY_FILE = 'assembly_graph_with_scaffolds.gfa'
CONTIG_FILE = 'contigs.fasta'

def str2bool(v):
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
		For example if tasks =[1 4] return [1 2 3 4]
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
		print("ERROR: There are more than two numbers in the task list!\n" + task_error_message)
		sys.exit()

	#valid_task_values = set(item.value for item in Pipeline_tasks)
	valid_task_values = [item.value for item in Pipeline_tasks]
	for task in tasks:
		if int(task) not in valid_task_values:
			print("ERROR: invalid task number(s)!\n" + task_error_message)
			sys.exit()

	if len(tasks)==2 and int(tasks[0])>int(tasks[1]):
		print("ERROR: The first task number should be smaller than the second task\
		 in the list!\n" + task_error_message)
		sys.exit()

	if len(tasks)==1 and int(tasks[0])==Pipeline_tasks.all.value:
		return valid_task_values
	if len(tasks)==1:
		return [int(tasks[0])]
	for task in list(range(int(tasks[0]), int(tasks[1])+1)):
		task_list.append(task)
	return task_list

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
	print("ERROR: "+message)
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
			if line.startswith('>'):
				continue
			return line

def insert_amr_in_location(genome_file, amr_seq, location, output_dir):
	"""
	To insert The AMR sequence in a given location (after a givel line) of the genome sequence
	Parameters:
		genome_file: the file in which the AMR file is to be inserted
		amr_seq:	the AMR sequence to be inserted
		location:	the line number to insert the AMR gene after it
		output_dir:	the output directory to store the result
	Return:
		the address of the generated file after the process
	"""
	#from /media/deta/salmonella.fasta to salmonella_1357.fasta when location=1357
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

def create_metagenome(ref_genome_files, number_of_insertions,insertion_type,
						insertion_locations, amr_file, output_dir):
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
	Return:
		The list of name of genome files with AMR inserted in them
	"""
	#some validations
	if len(params.number_of_insertions) != len(params.ref_genome_files):
		print("ERROR: Please specify the number of insertions for each reference genome.")
		sys.exit()
	if params.insertion_type == Insertion_type.assigned and \
		len(params.insertion_locations)!= sum(params.number_of_insertions):
		print("ERROR: Please specify the location of insertion for all insertions OR \
			if you prefer them to be chosen randomely choose 1 for --insertion_type")
		sys.exit()

	if not os.path.exists(output_dir):
		os.makedirs(output_dir)
	if (os.path.isfile(output_dir+'metagenome.fasta')):
		print('ERROR: A file named metagenome.fasta already exits in '+output_dir)
		sys.exit()

	insert_counter = 0
	genome_files = []
	amr_seq = retrieve_AMR(amr_file)
	for genome, insert_number in zip(ref_genome_files, number_of_insertions):
		if insertion_type.name == 'random':
			#find the number of lines in genome file
			lc = sum(1 for l in open(genome))
			#out = subprocess.Popen(['wc', '-l', genome], stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
			#stdout,stderr = out.communicate()
			#line_number = stdout.split()[0]
		for insertions in range(insert_number):
			if insertion_type.name == 'random':
				insert_location = random.randint(1,lc)
				new_file = insert_amr_in_location(genome, amr_seq, insert_location, output_dir)
			else:
				new_file = insert_amr_in_location(genome, amr_seq, insertion_locations[insert_counter], output_dir)
			genome_files.append(new_file)
			command = 'cat '+new_file+' >> ' + output_dir + 'metagenome.fasta'
			os.system(command)
			insert_counter+=1
	return genome_files

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
		print("ERROR: the read_length is not valid! In the current implementation it can\
		 be either 150 or 250!")
		sys.exit()

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
	os.system(command)

	#remove paths from GFA file
	gfa_file = output_dir + '/'+ ASSEMBLY_FILE
	command = "sed -i '/^P/d' " + gfa_file
	os.system(command)

	contig_file = output_dir + '/' + CONTIG_FILE

	return gfa_file, contig_file

def run_RGI(input_file, output_dir, seq_description, include_loose = False):
	"""
	"""
	rgi_dir = output_dir +"rgi_dir"
	if not os.path.exists(rgi_dir):
		os.makedirs(rgi_dir)
	output_file_name = rgi_dir +"/rgi_output_"+seq_description+"_"+datetime.datetime.now().strftime('%Y-%m-%d_%H-%M')
	command = "rgi main --input_sequence " + input_file + " --output_file " +\
		output_file_name +" --input_type protein --clean"
	if include_loose:
		command+=" --include_loose --exclude_nudge"
	os.system(command)
	seq_info_list = []
	with open(output_file_name + '.txt', newline = '') as rgi_file:
		rgi_reader = csv.reader(rgi_file, delimiter='\t')
		next(rgi_reader)
		for row in rgi_reader:
			seq_info = {'ORF_ID':row[0], 'gene':row[8].strip(),
			'prediction_type':row[5].strip(), 'best_identities':float(row[9])}
			# if row[2]=='+'
			# 	seq_info = {'gene':row[6].strip(),'length':row[2].strip(), 'start_pos':int(row[0]),
			# 	'end_pos':int(row[1]), 'prediction_type':row[3].strip(),'best_identities':float(row[7])}
			# elif row[2]=='-':
			# 	seq_info = {'gene':row[6].strip(),'length':row[2].strip(), 'start_pos':int(row[1]),
			# 	'end_pos':int(row[0]), 'prediction_type':row[3].strip(),'best_identities':float(row[7])}
			seq_info_list.append(seq_info)

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
	seq_file_name = output_dir+'tmp_'+datetime.datetime.now().strftime('%Y-%m-%d_%H-%M')+'.fasta'
	if os.path.isfile(seq_file_name):
		os.remove(seq_file_name)
	seq_file = open(seq_file_name, 'w')
	seq_file.write("> sequence:\n")
	seq_file.write(seq)
	seq_file.close()

	prokka_dir = output_dir+'prokka_dir_'+seq_description+'_'+datetime.datetime.now().strftime('%Y-%m-%d_%H-%M')
	if os.path.exists(prokka_dir):
		try:
			shutil.rmtree(prokka_dir)
		except OSError as e:
			print("Error: %s - %s." % (e.filename, e.strerror))

	command = prokka_prefix + 'prokka --metagenome --outdir '+\
		prokka_dir+' --prefix mygenome '+seq_file_name
	os.system(command)

	RGI_output_list = None
	if use_RGI:
		RGI_output_list = run_RGI(prokka_dir+"/mygenome.faa", output_dir, seq_description,
								RGI_include_loose)

	#analyze results and return something
	#key_values = ["gene", "length", "start_pos", "end_pos", "product"]
	#seq_info = dict.fromkeys(key_values)
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


	# #remove temporary files and folder
	# os.remove(seq_file_name)
	# try:
	# 	shutil.rmtree(prokka_dir)
	# except OSError as e:
	# 	print("Error: %s - %s." % (e.filename, e.strerror))

	return seq_info_list


def neighborhood_comparison(amr_file, ref_genome_files, neighborhood_seq_file, \
								contig_file, seq_length, output_dir, prokka_prefix,\
								use_RGI = True, RGI_include_loose = False):
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
		use_RGI:	RGI annotations incorporated for AMR annotation
		prokka_prefix: it's used to run prokka properly via docker or conda
	Return:
		the address of files stroing annotation information (annotation_detail_name,
			gene_file_name, product_file_name, visual_annoation)
	"""
	annotate_dir = output_dir+"annotation_"+str(seq_length)
	if os.path.exists(annotate_dir):
		try:
			shutil.rmtree(annotate_dir)
		except OSError as e:
			print("Error: %s - %s." % (e.filename, e.strerror))
	os.makedirs(annotate_dir)
	annotation_detail_name = annotate_dir+'/annotation_detail.csv'
	annotation_detail = open(annotation_detail_name, mode='w', newline='')
	#annotation_writer = csv.writer(annotation_detail, delimiter=',', quoting=csv.QUOTE_MINIMAL)
	annotation_writer = csv.writer(annotation_detail)
	if use_RGI:
		annotation_writer.writerow(["seq_name", "seq_value", "seq_length", "AMR_start_position",\
								"matched_contig", "gene", "prokka_gene_name", "product", \
								"length", "start_pos", "end_pos", 'RGI_prediction_type'])
	else:
		annotation_writer.writerow(["seq_name", "seq_value", "seq_length", "AMR_start_position",\
								"matched_contig", "gene", "product", "length", "start_pos", "end_pos"])
	gene_file_name = annotate_dir+'/seq_comparison_genes.txt'
	gene_file = open(gene_file_name, 'w')
	product_file_name = annotate_dir+'/seq_comparison_products.txt'
	product_file = open(product_file_name, 'w')

	#find neighborhood sequence in ref genome(s) and annotate using prokka
	amr_seq = retrieve_AMR(amr_file)
	ref_file_name = annotate_dir+'/ref_neighborhood_sequences_'+str(seq_length)+'_'+\
		datetime.datetime.now().strftime('%Y-%m-%d_%H-%M')+'.txt'
	ref_file = open(ref_file_name, 'w')
	for genome_file in ref_genome_files:
		genome_name = os.path.splitext(os.path.basename(genome_file))[0]
		ref_file.write("> " + genome_name + " :\n")
		seq, amr_start = extract_amr_neighborhood_in_ref_genome(amr_seq[:-1].lower(), genome_file, seq_length)
		ref_file.write(seq+"\n")
		seq_description = 'ref_'+genome_name
		seq_info_list = annotate_sequence(seq+"\n", seq_description, annotate_dir+'/',
										prokka_prefix, use_RGI, RGI_include_loose)
		myLine1 = myLine2 = seq_description + ':\t'
		for seq_info in seq_info_list:
			if use_RGI:
				annotation_writer.writerow([seq_description, seq, len(seq), amr_start, "", seq_info['gene'],
										seq_info['prokka_gene_name'], seq_info['product'],
										seq_info['length'],seq_info['start_pos'], seq_info['end_pos'],
										seq_info['RGI_prediction_type']])
			else:
				annotation_writer.writerow([seq_description, seq, len(seq), amr_start, "",
										seq_info['gene'], seq_info['product'],
										seq_info['length'],seq_info['start_pos'],
										seq_info['end_pos']])
			if seq_info['gene']=='':
				seq_info['gene']='UNKNOWN'
			myLine1+=seq_info['gene']+'---'
			myLine2+=seq_info['product']+'---'
		gene_file.write(myLine1[:-3]+'\n')
		product_file.write(myLine2[:-3]+'\n')
	ref_file.close()
	print("NOTE: The list of neighborhood sequences in reference genome(s) has\
	 		been stroed in " + ref_file_name)

	#annotate the sequences extraced from assembly graph
	counter = 1
	with open(neighborhood_seq_file, 'r') as read_obj:
		for line in read_obj:
			if line.startswith('>') or line.startswith('Path') or line.startswith('The'):
				continue
			#find out if there is a match for it in contig file
			contig_list = find_sequence_match(line, contig_file, annotate_dir, '_'+str(counter))
			contig_name = ""
			for contig in contig_list:
				if int(float(contig.identity))==100 and int(contig.matched_length) == (len(line)-1):
					contig_name = contig.name
			seq_description = 'extracted'+str(counter)
			seq_info_list = annotate_sequence(line, seq_description, annotate_dir+'/',
											prokka_prefix, use_RGI, RGI_include_loose)
			myLine1 = myLine2 = seq_description +':\t'
			amr_start = line.lower().find(amr_seq[:-1].lower())
			for seq_info in seq_info_list:
				if use_RGI:
					annotation_writer.writerow([seq_description, line[:-1], len(line[:-1]), amr_start+1,
										contig_name, seq_info['gene'], seq_info['prokka_gene_name'],
										seq_info['product'], seq_info['length'],
										seq_info['start_pos'], seq_info['end_pos'],
										seq_info['RGI_prediction_type']])
				else:
					annotation_writer.writerow([seq_description, line[:-1], len(line[:-1]),
										amr_start+1, contig_name, seq_info['gene'],
										seq_info['product'], seq_info['length'],
										seq_info['start_pos'], seq_info['end_pos']])
				if seq_info['gene']=='':
					seq_info['gene']='UNKNOWN'
				myLine1+=seq_info['gene']+'---'
				myLine2+=seq_info['product']+'---'
			gene_file.write(myLine1[:-3]+'\n')
			product_file.write(myLine2[:-3]+'\n')
			counter+=1

	annotation_detail.close()
	gene_file.close()
	product_file.close()
	print("NOTE: The comparison of neighborhood sequences are available in " +\
	 		annotation_detail_name+", "+gene_file_name+", "+product_file_name)

	visual_annoation = annotate_dir+'/gene_comparison.png'
	visualize_annotation(annotation_detail_name, output=visual_annoation)

	return annotation_detail_name, gene_file_name, product_file_name, visual_annoation


def main(params):
	#import pdb; pdb.set_trace()
	#Validate task values
	task_list = validate_task_values(params.task)

	graph_file =""
	seq_file= ""
	contigs_file = ""
	read1 = ""
	read2 = ""
	genome_amr_files = []
	gfa_file = None
	if Pipeline_tasks.metagenome_creation.value in task_list:
		genome_amr_files = create_metagenome(params.ref_genome_files, params.number_of_insertions,
					params.insertion_type, params.insertion_locations, params.amr_file,
					params.output_dir)

	if Pipeline_tasks.read_simulation.value in task_list:
		metagenome_file = params.output_dir + 'metagenome.fasta'
		read1, read2 = simulate_reads(metagenome_file, params.read_length, params.ART_PATH)

	if Pipeline_tasks.assembly.value in task_list:
		read1_file = verify_file_existence(read1, params.read1, \
			'please provide the address of both paired end reads')
		read2_file = verify_file_existence(read2, params.read2, \
			'please provide the address of both paired end reads')
		spade_output = params.output_dir + params.spades_output_dir
		graph_file, contigs_file = do_assembly(read1_file, read2_file, params.SPADES_PATH,
										spade_output, params.spades_thread_num,
										params.spades_error_correction)

	if Pipeline_tasks.graph_neighborhood.value in task_list:
		gfa_file = verify_file_existence(graph_file, params.gfa_file, \
				'please provide the address of the file containing the assembly graph')
		subgarph_file = neighborhood_graph_extraction(params.amr_file, gfa_file, params.graph_distance,
	 								params.output_dir, params.BANDAGE_PATH)
	if Pipeline_tasks.sequence_neighborhood.value in task_list:
		if not gfa_file:
			gfa_file = verify_file_existence(graph_file, params.gfa_file, \
					'please provide the address of the file containing the assembly graph')
		seq_file = neighborhood_sequence_extraction(params.amr_file, gfa_file, params.seq_length,
									params.output_dir, params.BANDAGE_PATH)

	if Pipeline_tasks.neighborhood_comparison.value in task_list:
		neighborhood_file = verify_file_existence(seq_file, params.neighborhood_seq_file,
				'please provide the address of the file containing all extracted\
					sequences from AMR neighborhood in the assembly graph')
		contig_file = verify_file_existence(contigs_file, params.contig_file,
				'please provide the address of the file containing contigs after assembly')
		if genome_amr_files:
			genome_files = genome_amr_files
		elif not params.genome_amr_files:
			print('ERROR: please provide the address of the files containing genome after AMR insertion')
			sys.exit()
		else:
			genome_files = params.genome_amr_files
		annotation_detail, gene_annotation, product_annotation, visual_annoation=\
			neighborhood_comparison(params.amr_file, genome_files,\
				neighborhood_file, contig_file, params.seq_length, params.output_dir,
				params.PROKKA_COMMAND_PREFIX, params.use_RGI, params.RGI_include_loose)


if __name__=="__main__":

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
	parser.add_argument('--amr_file','-A', type=str, default = params.amr_file,
		help = 'the path of the file containing the AMR gene sequence')
	#parser.add_argument('--ref_genome_number', type=int, default=params.ref_genome_number,
	#	help = 'the number of reference genomes that AMR genome will be inserted in them')
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
	parser.add_argument('--neighborhood_seq_file', type = str, default = params.neighborhood_seq_file,
		help = 'the address of the file containing all extracted neighborhood sequences in assembly graph')
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


	args = parser.parse_args()

	#updating params values
	params.task = args.task
	params.amr_file = args.amr_file
	#params.ref_genome_number = args.ref_genome_number
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
	params.neighborhood_seq_file = args.neighborhood_seq_file
	params.gfa_file = args.gfa_file
	params.contig_file = args.contig_file
	params.genome_amr_files = args.genome_amr_files
	params.read1 = args.read1
	params.read2 = args.read2
	params.spades_error_correction = args.spades_error_correction
	params.use_RGI = args.use_RGI
	params.RGI_include_loose = args.RGI_include_loose

	main(params)
