
"""
File:		amr_neighborhood_in_contigs.py
Aythor:		Somayeh Kafaie
Date:		March 2021
Purpose:	To find the neighborhood of AMRs in a contig file, compare them with
			that of the ref genomes and calculate the sentivity and precision

To run:
$ conda activate rgi
$ python amr_neighborhood_in_contigs.py
NOTE: It reads required parameters from params.py and the most important parameters need
to be set correctly there are:
params.seq_length, params.contig_file, params.amr_identity_threshold, params.amr_files,
params.ref_ng_annotations_file, params.main_dir, params.output_dir,
params.PROKKA_COMMAND_PREFIX, params.use_RGI,params.RGI_include_loose,
params.ref_genomes_available
NOTE: The result are available in the following directory:
params.output_dir+'contigs_output_'+str(params.seq_length)

"""

################################################################################

import os
import re
import argparse
import datetime
import csv
import subprocess
from csv import DictReader
import logging
import pandas as pd
import yaml

from sarand import utils, extract_neighborhood, find_amrs_in_sample
from sarand.utils import extract_files, retrieve_AMR, read_ref_annotations_from_db,\
		extract_up_down_from_csv_file, seqs_annotation_are_identical,\
		similar_seq_annotation_already_exist, split_up_down_info, annotate_sequence,\
		retreive_original_amr_name, extract_name_from_file_name, initialize_logger,\
		restricted_amr_name_from_modified_name, str2bool, validate_print_parameters_tools,\
		first_fully_covered_by_second
from sarand.extract_neighborhood import extract_amr_neighborhood_in_ref_genome
from sarand.find_amrs_in_sample import find_all_amrs_and_neighborhood

NOT_FOUND_FILE = 'not_found_amrs_in_contigs.txt'

# def read_ref_neighborhood_annotation(annotation_file):
# 	"""
# 	"""
# 	ref_up_info_list = []
# 	ref_down_info_list = []
# 	ref_amr_info_list = []
# 	with open(annotation_file, 'r') as myfile:
# 		myreader = DictReader(myfile)
# 		old_seq = ''
# 		seq_info =[]
# 		for row in myreader:
# 			if not row['seq_name'].startswith('extracted'):
# 				gene_info = {'seq_name':row['seq_name'], 'seq_value':row['seq_value'],
# 				 			'gene':row['gene'], 'length':row['length'],
# 							'start_pos':int(row['start_pos']),'end_pos':int(row['end_pos']),
# 							'target_amr':row['target_amr']}
# 				cur_seq = row['seq_name']
# 				if cur_seq!=old_seq:
# 					if (seq_info):
# 						amr_found, up_info, down_info, amr_info = extract_up_down_from_csv_file(seq_info)
# 						if amr_found:
# 							ref_amr_info_list.append(amr_info)
# 							if up_info and not similar_seq_annotation_already_exist(up_info, ref_up_info_list):
# 								ref_up_info_list.append(up_info)
# 							if down_info and not similar_seq_annotation_already_exist(down_info, ref_down_info_list):
# 								ref_down_info_list.append(down_info)
# 					seq_info = []
# 					old_seq = cur_seq
# 				seq_info.append(gene_info)
# 		amr_found, up_info, down_info, amr_info = extract_up_down_from_csv_file(seq_info)
# 		if amr_found:
# 			ref_amr_info_list.append(amr_info)
# 			if up_info and not similar_seq_annotation_already_exist(up_info, ref_up_info_list):
# 				ref_up_info_list.append(up_info)
# 			if down_info and not similar_seq_annotation_already_exist(down_info, ref_down_info_list):
# 				ref_down_info_list.append(down_info)
# 	return ref_amr_info_list, ref_up_info_list, ref_down_info_list

def extract_seq_neighborhood_and_annotate(amr_seq, amr_name, seq_length, contig_file,
								amr_threshold, contig_dir, prokka_prefix, use_RGI,
								RGI_include_loose, annotation_file_name):
	"""
	AMR alignment --> neighborhood extraction --> annotation for CAMI gold standard contig list
	"""
	up_info_list = []
	down_info_list = []
	amr_info_list = []
	seq_info_list = []
	restricted_amr_name = restricted_amr_name_from_modified_name(amr_name)
	seq_file_name = os.path.join(contig_dir, 'contig_sequences_'+restricted_amr_name+'_'+\
		str(seq_length)+'_'+datetime.datetime.now().strftime('%Y-%m-%d_%H-%M')+'.txt')
	seq_file = open(seq_file_name, 'w')
	#import pdb; pdb.set_trace()
	seq_list, contig_name_list = extract_amr_neighborhood_in_ref_genome(amr_seq, contig_file, seq_length, amr_threshold)
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
		seq_file.write("> " + contig_name+" :\n")
		seq_file.write(seq+"\n")
		seq_description = 'contig_'+contig_name
		seq_info = annotate_sequence(seq+"\n", seq_description, contig_dir,
									prokka_prefix, use_RGI, RGI_include_loose)
		myLine1 = myLine2 = seq_description + ':\t'
		found, amr_info , up_info, down_info, seq_info = split_up_down_info(seq, seq_info)
		if found:
			if up_info and not similar_seq_annotation_already_exist(up_info, up_info_list, contig_dir):
				up_info_list.append(up_info)
			if down_info and not similar_seq_annotation_already_exist(down_info, down_info_list, contig_dir):
				down_info_list.append(down_info)
			amr_info_list.append(amr_info)
		else:
			logging.error("no target amr was found in this contig sequence")
			#import pdb; pdb.set_trace()
			return [], [], [], []
		with open(annotation_file_name, 'a') as fd:
			writer = csv.writer(fd)
			for gene_info in seq_info:
				writer.writerow([seq_description, gene_info['seq_value'],
									len(gene_info['seq_value']),
									gene_info['gene'], gene_info['prokka_gene_name'],
									gene_info['product'], gene_info['length'],
									gene_info['start_pos'], gene_info['end_pos'],
									gene_info['RGI_prediction_type'],
									gene_info['family'], gene_info['target_amr']])
		seq_info_list.append(seq_info)
	seq_file.close()
	if not seq_list:
		logging.info(amr_name+' was not found in the contig list!')
	else:
		logging.info("NOTE: The list of neighborhood sequences in contigs has\
	 		been stroed in " + seq_file_name)
	return up_info_list, down_info_list, amr_info_list, seq_info_list

def evaluate_sequences_up_down(amr_name, summary_file, up_info_list, down_info_list,
								amr_info_list, ref_up_info_list, ref_down_info_list,
								ref_amr_info_list, found_file):
	"""
	To compare upstream/downstream annotations from contig with those of ref genomes
	"""
	out_dir = os.path.dirname(summary_file)
	ref_len =  len(ref_up_info_list)+len(ref_down_info_list)
	#if amr was not found in the ref genomes
	if ref_len==0 and len(ref_amr_info_list)==0:
		with open(summary_file,'a') as fd:
			writer = csv.writer(fd)
			#writer.writerow([amr_name, 0, 0, 0, len(up_info_list)+len(down_info_list),-1, -1])
			writer.writerow([amr_name, 0, 0, 0, len(up_info_list)+len(down_info_list),-1, -1, True])
		logging.error(amr_name+" was not found in the ref genomes!!!")
		import pdb; pdb.set_trace()
		return -1, -1
	#If AMR was not found in the contig list
	if not amr_info_list and not up_info_list and not down_info_list:
		with open(found_file, 'a') as fa:
			fa.write(amr_name+'\n')
		with open(summary_file,'a') as fd:
			writer = csv.writer(fd)
			writer.writerow([amr_name, 0, 0, ref_len, 0,0, 0, True])
			#writer.writerow([amr_name, 0, 0, ref_len, 0,0, 0])
		return 0, 0
	#find the number of unique true-positives, all false positives, total found cases, all unique true cases
	unique_tp = 0
	for ref_info in ref_up_info_list:
		for seq_info in up_info_list:
			#if seqs_annotation_are_identical(ref_info, seq_info, out_dir):
			#if all ref genes are available in the extracted one in the same order
			# and it doesn't matter if extracted one has more genes!
			if first_fully_covered_by_second(ref_info, seq_info, out_dir,
										in_reversed_order = True):
				unique_tp+=1
				break
	for ref_info in ref_down_info_list:
		for seq_info in down_info_list:
			#if seqs_annotation_are_identical(ref_info, seq_info, out_dir):
			#if all ref genes are available in the extracted one in the same order
			# and it doesn't matter if extracted one has more genes!
			if first_fully_covered_by_second(ref_info, seq_info, out_dir):
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
				if seqs_annotation_are_identical(ref_info, seq_info, out_dir):
					found_identical_annotation = True
					break
			if not found_identical_annotation:
				false_positive+=1
	for seq_info in down_info_list:
		if seq_info:
			found_cases_len+=1
			found_identical_annotation = False
			for ref_info in ref_down_info_list:
				if seqs_annotation_are_identical(ref_info, seq_info, out_dir):
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
		# writer.writerow([amr_name, unique_tp, false_positive, ref_len, found_cases_len,
		# 				sensitivity, precision])
		writer.writerow([amr_name, unique_tp, false_positive, ref_len, found_cases_len,
						sensitivity, precision, False])

	return sensitivity, precision

def annotate_sequence_bundle(contig_ng_info, out_dir, prokka_prefix, use_RGI,
								RGI_include_loose):
	"""
	"""
	annotate_dir = os.path.join(out_dir,'contig_annotations')
	if not os.path.exists(annotate_dir):
		os.makedirs(annotate_dir)
	for amr_info_list in contig_ng_info:
		amr_name = amr_info_list[0]
		restricted_amr_name = restricted_amr_name_from_modified_name(amr_name)
		# mydir = annotate_dir+restricted_amr_name
		# if not os.path.exists(mydir):
		# 	os.makedirs(mydir)
		annotation_file_name =os.path.join(annotate_dir, 'contig_annotation_'+restricted_amr_name+'.csv')
		with open(annotation_file_name, 'a') as fd:
			writer = csv.writer(fd)
			writer.writerow(['seq_name', 'seq_value', 'seq_length', 'gene', 'prokka_gene_name',
							'product', 'length', 'start_pos', 'end_pos', 'RGI_prediction_type',
							 'family', 'target_amr'])
		for index, amr_info in enumerate(amr_info_list[1]):
			contig_name = amr_info['contig']
			seq = amr_info['seq']
			seq_description = 'contig_'+amr_name+'_'+contig_name.replace(' ','_').replace('.','').replace(',','')
			annotation_prefix = 'contig_'+ restricted_amr_name +'__'+str(index)
			seq_info = annotate_sequence(seq+"\n", annotation_prefix, annotate_dir,
                    prokka_prefix, use_RGI, RGI_include_loose)
			found, _ , _, _, _ = split_up_down_info(seq, seq_info)
			if not found:
				logging.error("no target amr was found in this contig sequence: "+contig_name)
			with open(annotation_file_name, 'a') as fd:
				writer = csv.writer(fd)
				for gene_info in seq_info:
					writer.writerow([seq_description, gene_info['seq_value'],
										len(gene_info['seq_value']),
										gene_info['gene'], gene_info['prokka_gene_name'],
										gene_info['product'], gene_info['length'],
										gene_info['start_pos'], gene_info['end_pos'],
										gene_info['RGI_prediction_type'],
										gene_info['family'], gene_info['target_amr']])
		logging.info("NOTE: The annotation of neighborhood sequences in contigs for "+\
			amr_name+"has been stroed in " + annotation_file_name)

def find_contig_amrs_main(params):
	"""
	"""
	#creating a directory for results
	contig_dir = os.path.join(params.output_dir, 'contigs_output_'+str(params.seq_length))
	if not os.path.exists(contig_dir):
		os.makedirs(contig_dir)
	found_file = os.path.join(contig_dir , NOT_FOUND_FILE)

	if params.ref_genomes_available:
		#extracting all amr files
		amr_files = extract_files(params.amr_files, 'please provide the address of the AMR gene(s)')
		# create a db from contig file
		db_command = subprocess.run(["makeblastdb","-in", params.contig_file, "-parse_seqids",
									"-dbtype", "nucl"], stdout=subprocess.PIPE, check= True)
		logging.info(db_command.stdout.decode('utf-8'))
		# command = 'makeblastdb -in '+params.contig_file +' -parse_seqids -dbtype nucl'
		# os.system(command)
		#to find the annotation of ref genomes for all AMRs
		#set up the file to store the summary metrics
		summary_file = os.path.join(contig_dir, 'contig_summaryMetrics_up_down_'+\
		datetime.datetime.now().strftime('%Y-%m-%d_%H-%M')+'.csv')
		with open(summary_file,'a') as fd:
			writer = csv.writer(fd)
			writer.writerow(['AMR', 'Unique_TP#', 'FP#', 'Unique_True#', 'found#','sensitivity', 'precision', 'not_found'])
		df = pd.read_csv(params.ref_ng_annotations_file, skipinitialspace=True,  keep_default_na=False)
		amr_groups = df.groupby('target_amr')
		#Going over amrs one by one to extract their neighborhood and annotate and evaluate
		average_precision = 0
		average_sensitivity = 0
		for amr_file in amr_files:
			restricted_amr_name = extract_name_from_file_name(amr_file)
			mydir = os.path.join(contig_dir, restricted_amr_name)
			amr_seq, amr_name = retrieve_AMR(amr_file)
			if not os.path.exists(mydir):
				os.makedirs(mydir)
			annotation_file_name =os.path.join(mydir, 'contig_annotation_'+restricted_amr_name+'_'+str(params.seq_length)+'.csv')
			with open(annotation_file_name, 'a') as fd:
				writer = csv.writer(fd)
				writer.writerow(['seq_name', 'seq_value', 'seq_length', 'gene', 'prokka_gene_name',
								'product', 'length', 'start_pos', 'end_pos', 'RGI_prediction_type',
								 'family', 'target_amr'])
			#extract neighborhoods in contigs and annotate
			up_info_list, down_info_list, amr_info_list, seq_info_list =\
				extract_seq_neighborhood_and_annotate(amr_seq, amr_name, params.seq_length,
									params.contig_file, params.amr_identity_threshold, mydir,
									params.PROKKA_COMMAND_PREFIX, params.use_RGI,
									params.RGI_include_loose, annotation_file_name)
			#Read ref neighborhood annotations from corresponding file
			ref_up_info_list, ref_amr_info_list, ref_down_info_list =\
				read_ref_annotations_from_db(amr_groups, amr_name)
			original_amr_name = retreive_original_amr_name(amr_name)
			sensitivity, precision = evaluate_sequences_up_down(original_amr_name, summary_file,
										up_info_list, down_info_list, amr_info_list,
										ref_up_info_list, ref_down_info_list, ref_amr_info_list,
										found_file)
			average_precision+=precision
			average_sensitivity+=sensitivity
			logging.info('For "'+amr_name+'": sensitivity= '+str(sensitivity)+' precision = '+ str(precision))
		logging.info('average precision: '+str(average_precision/len(amr_files))+'\naverage sensitivity: '+ str(average_sensitivity/len(amr_files)))
	else:
		#read card sequences and do neighborhood extraction and annotations for blast hits
		contig_ng_info = find_all_amrs_and_neighborhood(params.amr_db, params.contig_file,
										contig_dir, params.seq_length,
										params.amr_identity_threshold, type='contig')
		annotate_sequence_bundle(contig_ng_info, contig_dir, params.PROKKA_COMMAND_PREFIX,
									params.use_RGI, params.RGI_include_loose)
	logging.info("neighborhood extraction and annotation from cotigs is done!")

def update_contig_params(params, config):
	"""
	"""
	config_params = config.keys()
	main_dir_changed = False
	if 'main_dir' in config_params and os.path.realpath(config['main_dir'])!=os.path.realpath(params.main_dir):
		main_dir_changed = True
		#changes params variables dependant on main_dir only accessible through params.py
		params.output_dir = params.output_dir.replace(params.main_dir.rstrip(' /'), config['main_dir'].rstrip(' /'))
		params.amr_files = params.amr_files.replace(params.main_dir.rstrip(' /'), config['main_dir'].rstrip(' /'))
		params.ref_ng_annotations_file = params.ref_ng_annotations_file.replace(params.main_dir.rstrip(' /'), config['main_dir'].rstrip(' /'))
		params.metagenome_file = params.metagenome_file.replace(params.main_dir.rstrip(' /'), config['main_dir'].rstrip(' /'))
		params.ng_seq_files = params.ng_seq_files.replace(params.main_dir.rstrip(' /'), config['main_dir'].rstrip(' /'))
		params.ng_path_info_files = params.ng_path_info_files.replace(params.main_dir.rstrip(' /'), config['main_dir'].rstrip(' /'))
	if 'PROKKA_COMMAND_PREFIX' in config_params:
		params.PROKKA_COMMAND_PREFIX = config['PROKKA_COMMAND_PREFIX']
	if 'main_dir' in config_params:
		params.main_dir = config['main_dir']
	if 'amr_db' in config_params and config['amr_db']!='':
		params.amr_db = config['amr_db']
	if 'seq_length' in config_params:
		params.seq_length = config['seq_length']
	if 'contig_file' in config_params:
		params.contig_file = config['contig_file']
	elif main_dir_changed:
		params.contig_file =params.contig_file.replace(params.main_dir.rstrip(' /'), config['main_dir'].rstrip(' /'))
	if 'amr_identity_threshold' in config_params:
		params.amr_identity_threshold = config['amr_identity_threshold']
	if 'ref_genomes_available' in config_params:
		params.ref_genomes_available = config['ref_genomes_available']
	if 'use_RGI' in config_params:
		params.use_RGI = config['use_RGI']
	if 'RGI_include_loose' in config_params:
		params.RGI_include_loose = config['RGI_include_loose']
	if main_dir_changed:
		params.main_dir = config['main_dir']

	return params

if __name__=="__main__":
	import params
	text = 'This code is used to find the AMR neighborhood in the contig file'
	parser = argparse.ArgumentParser(description=text)
	parser.add_argument('--config_file', '-C', type = str, default='',
		help = 'the config file to set parameters for find_contig_amrs()')
	args = parser.parse_args()
    # Read config file into a dictionery
	print("Reading the config file '"+args.contig_file+"' ...")
	with open(args.config_file, 'r') as yamlfile:
		data = yaml.load(yamlfile, Loader=yaml.FullLoader)
	params = update_contig_params(params, data)
	log_name = 'logger_'+datetime.datetime.now().strftime('%Y-%m-%d_%H-%M')+'.log'
	initialize_logger(args.main_dir, log_name)
	validate_print_parameters_tools(params, "find_contig_amrs")
	find_contig_amrs_main(args)
