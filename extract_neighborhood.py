"""
File:		extract_neighborhood.py
Aythor:		Somayeh Kafaie
Date:		July 2020
Purpose:	To extract the neighborhood of an AMR gene from an assembly graph

To run:
- The Bandage+BLAST implementation (type = 1):
	python extract_neighborhood.py --amr/-A <AMR gene file path in FASTA format>
	--gfa/-G <GFA assembly graph> --distance/-D <max distance from AMR gene (default=1)>
	--length/-L <length of the linear sequence around AMR gene to be extracted (default = 1000)>

"""

################################################################################

import sys
import os
import errno
import gfapy
import re
import argparse
import difflib
import datetime
import csv
import collections
from Bio import SeqIO
from gfapy.sequence import rc
import shutil
import logging
import enum

from utils import reverse_sign, find_node_orient, find_node_name, find_node_name_orient,\
					exist_in_path
from params import Assembler_name
#from find_seq_in_contigs import find_sequence_match

BANDAGE_PATH = '/media/Data/tools/Bandage_Ubuntu_dynamic_v0_8_1/Bandage'
OUT_DIR = 'output/'
BOTH_DIR_RECURSIVE = True
#MAX_K_MER_SIZE = 55
#if this feature is true after that the extracted length reaches a threshold, we don't go
#over all edges in the neighbourhood (in pre_sequence and post_sequence extraction)
#and only extract the sequence from the longest one to reduce the processing time!!!
ONLY_LONGEST_EDGE = False
TEMP_DIR = 'temp/'

def retrieve_AMR(file_path):
	"""
	To read the AMR gene from the fasta file.
	Parameters:
		file_path:	the address of the file containing the AMR gene
	Return:
		the sequence of the AMR gene in lower case
	"""
	with open(file_path) as fp:
		for line in fp:
			if line.startswith('>'):
				continue
			return line.lower()

def compare_two_sequences(seq1, seq2, output_dir, threshold = 90, blast_ext = ''):
	"""
	To compare one sequence (shorter sequence) against the other one (longer sequence) using blastn
	"""
	#make sure seq1 is the longer sequence
	if len(seq1)<len(seq2):
		seq1, seq2 = seq2, seq1
	#write the query sequence into a fasta file
	query_file_name = output_dir+'query.fasta'
	with open(query_file_name, 'w') as query_file:
		query_file.write('> query \n')
		query_file.write(seq2)
	#write the query sequence into a fasta file
	subject_file_name = output_dir+'subject.fasta'
	with open(subject_file_name, 'w') as subject_file:
		subject_file.write('> subject \n')
		subject_file.write(seq1)
	#run blast query for alignement
	blast_file_name = output_dir+'blast'+blast_ext+'.csv'
	command = 'blastn -query '+query_file_name+' -subject '+subject_file_name+\
		' -outfmt 10  > '+ blast_file_name
	os.system(command)

	with open(blast_file_name, 'r') as file1:
		myfile = csv.reader(file1)
		for row in myfile:
			identity=int(float(row[2]))
			coverage = int(float(row[3])/len(seq2)*100)
			if identity>=threshold and coverage>=threshold:
				return True

	return False

def similar_sequence_exits(seq_list, query_seq, output_dir):
	"""
	To check if any similar sequence to query_seq exists in seq_list
	if it finds any similar one with a shorter length, replaces that with query_Seq
	Parameters:
		seq_list: 	the list of sequences
		query_seq:	the query sequence
		output_dir:	the output directory
	Return:
		the index of the similar seq if it has a shorter length than query_seq (to be replace)
		and
		True if a seq similar to query_seq was found in seq_list

	"""
	for i, seq in enumerate(seq_list):
		if compare_two_sequences(seq, query_seq, output_dir, blast_ext = str(i)):
			#if the new sequence has a longer length replace the shorter similar oe with it
			if len(seq)<len(query_seq):
				return i, True
			else:
				return -1, True
	return -1, False

def find_sequence_match(query, contig_file, out_dir = TEMP_DIR, blast_ext = '',
						threshold =95, max_target_seqs = 100):
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
	# #Create DB from contig file
	# command = 'makeblastdb -in '+contig_file +' -parse_seqids -dbtype nucl'
	# os.system(command)
	#write the sequence into a fasta file
	query_file = out_dir+'/query.fasta'
	file = open(query_file, 'w')
	file.write(query)
	file.close()
	#run blast query for alignement
	blast_file_name = out_dir+'/blast'+blast_ext+'.csv'
	command = 'blastn -query '+query_file+' -db '+contig_file+\
		' -task blastn -outfmt 10 -max_target_seqs '+str(max_target_seqs)+' -evalue 0.5 -perc_identity '+str(threshold-1)+' > '+ blast_file_name
	# command = 'blastn -query '+query_file+' -db '+contig_file+\
	# 	' -task blastn -outfmt 10 -max_target_seqs 5 -evalue 0.5 -perc_identity 95 > '+ blast_file_name

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

def extract_amr_neighborhood_in_ref_genome(amr_seq, ref_path, neighborhood_len, threshold = 95):
	"""
	To find the neighborhood of amr_seq in the genome (ref_path)
	Parameters:
		amr_seq:	the sequence of the AMR gene
		ref_path:	the address of the reference genome
		neighborhood_len: 	the length of neighborhood from each side to be extracted
		threshold: used for identity and coverage
	Return:
		the neighborhood(s) of amr_seq with maximum length = 2*len + length(amr_seq)
			and the start position of the AMR in the extracted sequence(s)
	# In case AMR exists in multiple locations of the genome this function returns
		multiple neighborhood sequences
	"""
	# copy ref_path in a new temp location and folder that can be deleted at the end
	temp_dir = 'tmp_ref_path'
	if not os.path.exists(temp_dir):
		try:
			os.makedirs(temp_dir)
		except OSError as exc:
			if exc.errno != errno.EEXIST:
				raise
			pass
	# tmp_path = temp_dir + '/'+os.path.basename(ref_path)
	# shutil.copyfile(ref_path, tmp_path)

	# check if the sequence can be found in the genome
	#record_list = find_sequence_match(amr_seq, tmp_path, temp_dir, threshold = threshold)
	record_list = find_sequence_match(amr_seq, ref_path, temp_dir, threshold = threshold)

	# to get ride of \n
	len_amr_seq = len(amr_seq)-1

	# read record_list and store cases with identity and coverage above the threshold
	# and extract the beginning and end of AMR in the ref_file
	found = False
	amr_start_list = []
	amr_end_list = []
	record_name_list = []
	# find the start and end of AMR for all found cases above the threshold
	for record in record_list:
		if int(float(record.identity))>=threshold and int(float(record.matched_length)/len_amr_seq*100)>=threshold:
			found = True
			amr_start_list.append(int(record.c_start))
			amr_end_list.append(int(record.c_end))
			record_name_list.append(record.name)
			# amr_start = int(record.c_start)
			# amr_end = int(record.c_end)
			# record_name = record.name
			# if (amr_start not in amr_start_list) and (amr_end not in amr_end_list):
			# 	amr_start_list.append(amr_start)
			# 	amr_end_list.append(amr_end)
			# 	record_name_list.append(record_name)

	# extract neighborhood sequence(s)
	seq_list = []
	contig_name_list = []
	if found:
		for i, amr_start in enumerate(amr_start_list):
			amr_end = amr_end_list[i]
			record_name = record_name_list[i]
			seq = ''
			reverse_complement = False
			# extract sequence from both sides
			if amr_start > amr_end:
				amr_start, amr_end = amr_end, amr_start
				reverse_complement = True
			for record in SeqIO.parse(open(ref_path,'r'),'fasta'):
				#if record.id == record_name:
				if record_name == record.id:
					amr_start = amr_start -1
					if amr_start-neighborhood_len >= 0:
						#seq = str(record.seq[amr_start-neighborhood_len:amr_end])
						seq = str(record.seq[amr_start-neighborhood_len:amr_start]).upper()+str(record.seq[amr_start:amr_end]).lower()
						amr_start = neighborhood_len
					else:
						#seq = str(record.seq[:amr_end])
						seq = str(record.seq[:amr_start]).upper()+str(record.seq[amr_start:amr_end]).lower()
					if (amr_end-1+neighborhood_len)<=len(record):
						seq+=str(record.seq[amr_end:amr_end+neighborhood_len])
					else:
						seq+=str(record.seq[amr_end:])
					if reverse_complement:
						seq = rc(seq)
					seq_list.append(seq)
					contig_name_list.append(record_name)
					break

	# delete the temp folder
	if os.path.exists(temp_dir):
		try:
			shutil.rmtree(temp_dir)
		except OSError as e:
			logging.error("Error: %s - %s." % (e.filename, e.strerror))
			#print("Error: %s - %s." % (e.filename, e.strerror))
	#added by 1 because the index starts from 0 in the string
	return seq_list, contig_name_list

def extract_neighbourhood(segment, segment_list, edge_list, distance):
	"""
	To extract the neighbourghood of a node
	Parameters:
		segment: the node we are looking for its neighborhood
		segment_list: the list of segments in the neighbourhood so far
		edge_list: the list of edges in the neighbourhood so far
		distance: the maximum distance of neighbourhood we are going to extract from AMR gene
	Return:
		updated edge_list and segment_list
	"""
	if distance <= 0:
		return segment_list, edge_list
	for link in reversed(segment.dovetails_L):
		if link.from_segment.name != segment.name:
			other_side = link.from_segment
		elif link.to_segment.name != segment.name:
			other_side = link.to_segment
		else:
			logging.info("WARNING: the link is a loop for node#"+segment.name)
			#print("WARNING: the link is a loop for node#"+segment.name)
			other_side = segment
		if link not in edge_list:
			edge_list.append(link)
		segment_list.add(other_side.name)
		segment_list, edge_list =\
			extract_neighbourhood(other_side, segment_list, edge_list, distance - 1)
	for link in reversed(segment.dovetails_R):
		if link.from_segment.name != segment.name:
			other_side = link.from_segment
		elif link.to_segment.name != segment.name:
			other_side = link.to_segment
		else:
			logging.info("WARNING: the link is a loop for node#"+segment.name)
			#print("WARNING: the link is a loop for node#"+segment.name)
			other_side = segment
		if link not in edge_list:
			edge_list.append(link)
		segment_list.add(other_side.name)
		segment_list, edge_list =\
			extract_neighbourhood(other_side, segment_list, edge_list, distance - 1)
	return segment_list, edge_list

def create_graph(myGraph, segment_list, edge_list):
	"""
	To create a graph from the list of segments and edges
	Parameters:
		myGraph: The original graph (the new graph is a subgraph of this graph)
		segment_list: the list of segments of the new graph
		edge_list: the list of edges of the new graph
	Return:
		the created graph which is a subgraph of myGraph
	"""
	subGraph = gfapy.Gfa(version='gfa1')
	for segment_name in segment_list:
		seg = myGraph.segment(segment_name)
		seg.disconnect()
		subGraph.append(seg)
	for edge in edge_list:
		if edge.is_connected():
			edge.disconnect()
		subGraph.append(edge)
	return subGraph

def extract_subgraph(nodes, distance, file_path):
	"""
	To extract the neighborhood.
	This function calls the recursive function "extract_neighbourhood"
	for the first and last nodes in the path to do the job
	Parameters:
		nodes: the list of nodes we are interested in their neighbourhood
		distance: the maximum distance of neighbourhood we are going to extract from AMR gene
		file_path: the address of the original assembly graph in GFA format
	Return:
		the subgraph representing the neighbourhood of the node_lists in a given max distance
	"""
	myGraph = gfapy.Gfa.from_file(file_path)
	segment_list=set()
	edge_list = []
	#for now we just find the neighborhood of the first and last nodes in the path
	#Add first node in AMR path and find its neighborhood
	s = myGraph.segment(re.sub('[+-]', '', nodes[0]))
	segment_list.add(s.name)
	segment_list, edge_list = extract_neighbourhood(s, segment_list, edge_list, distance)
	#Add last node in AMR path and find its neighborhood
	if (len(nodes) > 1):
		segment_list_last = set()
		edge_list_last = []
		s = myGraph.segment(re.sub('[+-]', '', nodes[-1]))
		segment_list.add(s.name)
		segment_list_last, edge_list_last = extract_neighbourhood(s, segment_list_last, edge_list_last, distance)
		#merge the lists
		segment_list = segment_list.union(segment_list_last)
		for edge in edge_list_last:
			if edge not in edge_list:
				edge_list.append(edge)

	return myGraph, segment_list, edge_list

def extract_nodes_in_path(path):
	"""
	Parameters:
		path:	a list of nodes with -/+ tail and comma separated : e.g., '(1363) 69-, 2193+ (1786)'
	Return:
		node_list:	list of node numbers -> e.g., [69, 2193]
		orientation_list: list of orientation of nodes -> e.g., [-, +]
		start_pos:	where in the first node the sequence has started -> e.g., 1363
		end_pos:	where in the last node the sequence ended  -> e.g., 1786
	"""
	start_pos = 0
	end_pos = 0
	if path.startswith('('):
		index = path.find(')')
		start_pos = int(path[1:index])
	if path.endswith(')'):
		index = path.rfind('(')
		end_pos = int(path[index+1:-1])
	#Remove text between ()
	path = (re.sub("\((.*?)\)", "", path)).strip()
	node_list = []
	orientation_list = []
	nodes = path.split(',')
	for node in nodes:
		if '-' in node:
			orientation_list.append('-')
		else:
			orientation_list.append('+')
		node = re.sub('[+-]', '', node.split()[0])
		node_list.append(node)
	return node_list, orientation_list, start_pos, end_pos

def write_sequences_to_file(sequence_list, path_list, file_name):
	"""
	To write the extracted sequences and paths in neighborhood of the AMR gene in a file
	Parameters:
		sequence_list: the list of sequences extracted in AMR gene neighborhood
		path_list: the list of all paths (a list of nodes) in AMR gene neighborhood
		file_name: the name of file to be writtn in
	"""
	file = open(file_name, 'a+')
	for seq, path in zip(sequence_list, path_list):
		myLine = "Path-> " + ", ".join(path) + ":"
		file.write(myLine+"\n")
		file.write(seq+"\n")
	file.close()

def sequence_on_orientation(seq, orient):
	"""
	To	return the sequence based on the orientation used in the graph
	Parameters:
	 	orient:	+ / -
		seq:	the original sequence
	Return: reverse complement for -; otherwise seq itself
		for orient = '-' and seq = 'AACTG' -> 'CAGTT'
	"""
	if orient=='-':
			return rc(seq)
	else:
		return seq

def find_overlap(from_node, to_node, from_orient, to_orient, validate_overlap_size = True):
	"""
	To return the size of overlap beween  from_node and to_node
	e.g.; from_node = 2, to_node = 3, from_orient = +, to_orient= +, edge = [L 2 + 3 + 55M]
				---->  return 55
	Parameters:
		from_node:	the first node of the edge
		to_node:	the second node of the edge
		from_orient:the orientation of the first node of the edge
		to_orient:	the orientation of the second node of the edge
		validate_overlap_size: to verify if the extracted size for overlap is correct and matches
	Return:
		the size of overlap
	"""
	size = 0
	for edge in (from_node.dovetails_L + from_node.dovetails_R):
		# 3561+ --> 15082+
		if edge.from_segment == from_node and edge.to_segment == to_node and\
			edge.from_orient == from_orient and edge.to_orient == to_orient:
			if edge.overlap[0].code == 'M':
				size = edge.overlap[0].length
				if (size>0 and validate_overlap_size):
					seq1 = sequence_on_orientation(edge.from_segment.sequence, edge.from_orient)
					seq2 = sequence_on_orientation(edge.to_segment.sequence, edge.to_orient)
					if seq1[len(seq1)-size:] != seq2[:size]:
						logging.info("WARNING: the overlap doesn't match: "+str(edge))
						logging.info("first sequence: "+seq1)
						logging.info("second sequence: "+seq2)
						# print("WARNING: the overlap doesn't match: "+str(edge))
						# print("first sequence: "+seq1)
						# print("second sequence: "+seq2)
			else:
				logging.info("WARNING: the sequence from node "+from_node.name+" to node "+
						to_node.name + "has "+ edge.overlap[0] + " instead of overlap!")
				# print("WARNING: the sequence from node "+from_node.name+" to node "+
				# 		to_node.name + "has "+ edge.overlap[0] + " instead of overlap!")
			return size
		#sometimes we have the reverse complement of the edge and not itself
		# 15082- --> 3561-
		elif edge.from_segment == to_node and edge.to_segment == from_node and\
			edge.from_orient == reverse_sign(to_orient) and edge.to_orient == reverse_sign(from_orient):
			if edge.overlap[0].code == 'M':
				size = edge.overlap[0].length
				if (validate_overlap_size):
					seq1 = sequence_on_orientation(edge.to_segment.sequence, reverse_sign(edge.to_orient))
					seq2 = sequence_on_orientation(edge.from_segment.sequence, reverse_sign(edge.from_orient))
					if seq1[len(seq1)-size:] != seq2[:size]:
						logging.info("WARNING: the overlap doesn't match: "+str(edge))
						logging.info("first sequence: "+seq1)
						logging.info("second sequence: "+seq2)
						# print("WARNING: the overlap doesn't match: "+str(edge))
						# print("first sequence: "+seq1)
						# print("second sequence: "+seq2)
			else:
				logging.info("WARNING: the sequence from node "+from_node.name+" to node "+
						to_node.name + "has "+ edge.overlap[0] + " instead of overlap!")
				# print("WARNING: the sequence from node "+from_node.name+" to node "+
				# 		to_node.name + "has "+ edge.overlap[0] + " instead of overlap!")
			return size
	return size

def remove_identical_sub_sequences(sequence_list, path_list):
	"""
	To keep only one instance in case of multiple identical extracted sequences or
	when a sequence is a subsequence of another one
	Parameters:
		sequence_list:	list of extracted sequences
		path_list:		list of extracted paths
	Return sequence_list and path_list in which every sequence/path is unique and
		none of them is subssequence/subpath of another one
	"""
	#In case two sequences are identical keep the one with the lower index in the list
	identical_sub_set = set()
	for i, seq1 in enumerate(sequence_list):
		for j, seq2 in enumerate(sequence_list):
			if i < j:
				if seq1==seq2 or seq2 in seq1:
					identical_sub_set.add(j)
				elif seq1 in seq2:
					identical_sub_set.add(i)

	identical_sub_list = sorted(identical_sub_set)
	for index in reversed(identical_sub_list):
		del sequence_list[index]
		del path_list[index]

	return sequence_list, path_list

def reverse_path(path):
	"""
	To generate the reverse of a path;
		e.g., [1+, [8-, 12+], 9-] --> [9+, [12-, 8+,] 1-]
	"""
	mypath = []
	items = ['[', ']', '{', '}']
	reversed_items = [']', '[', '}', '{']
	for node in reversed(path):
		mynode = ''
		num=''
		for ch in node:
			if ch in items:
				mynode = reversed_items[items.index(ch)] + mynode
			elif ch=='-' or ch=='+':
				num+=reverse_sign(ch)
				mynode = num + mynode
			else:
				num+=ch
		mypath.append(mynode)
	return mypath

# def generate_all_paths(pre_path_list, post_path_list, node_dir_list, path_dir):
# 	"""
# 	Generate the list of all paths which are combination of the input lists
# 	e.g., pre_path_list = [['1-','2+']], post_path_list=[['7+'],['8+','9+','10-']],
# 	node_dir_list = ['[13+','14-]']
# 	- if path_dir = 'same' --> 		[
# 						['1-','2+','[13+','14-]','7+'],
# 						['1-','2+','[13+','14-]','8+','9+','10-']
# 									]
# 	- if path_dir = 'reverse' --> 	[
# 						['7-', '[14+', '13-]', '2-', '1+'],
# 						['10+', '9-', '8-', '[14+', '13-]', '2-', '1+']
# 									]
# 	Parameters:
# 		pre_path_list: the list of nodes (and their orientation) preceding the AMR node
# 		post_path_list: the list of nodes (and their orientation) following the AMR node
# 		node_dir_list: the list of nodes (and their orientation) presenting the AMR node
# 		path_dir: 'same' if the AMR path direction (node_dir_list) is the same as what
# 			blast returned (presenting actual AMR sequence) and 'reverse' if the AMR path
# 			direction (node_dir_list) is the reverse of blast result (presenting the
# 			reverse complement of AMR sequence); in the latter, we need to return the
# 			reverse of the path to have the actual AMR sequence and not its reverse.
# 	"""
# 	path_list = []
# 	if pre_path_list and post_path_list:
# 		for pre_path in pre_path_list:
# 			for post_path in post_path_list:
# 				if path_dir == 'same':
# 					path_list.append(pre_path +node_dir_list + post_path)
# 				elif path_dir == 'reverse':
# 					path_list.append(reverse_path(pre_path +node_dir_list + post_path))
# 	elif pre_path_list:
# 		for pre_path in pre_path_list:
# 			if path_dir == 'same':
# 				path_list.append(pre_path +node_dir_list)
# 			elif path_dir == 'reverse':
# 				path_list.append(reverse_path(pre_path +node_dir_list))
# 	elif post_path_list:
# 		for post_path in post_path_list:
# 			if path_dir == 'same':
# 				path_list.append(node_dir_list + post_path)
# 			elif path_dir == 'reverse':
# 				path_list.append(reverse_path(node_dir_list + post_path))
# 	else:
# 		if path_dir == 'same':
# 			path_list.append(node_dir_list)
# 		elif path_dir == 'reverse':
# 			path_list.append(reverse_path(node_dir_list))
# 	return path_list

def extract_found_amr(myGraph, node_list, orientation_list, start_pos, end_pos):
	"""
	To extract the sequence representing AMR gene from the path found by
	Bandage+BLAST presented in node_list; I could have used the original AMR gene
	sequence but I thought there might be cases that the assembler can't assemble
	the entire AMR gene or there are some errors in assembled one.
	Parameters:
		myGraph: The original graph
		node_list:	the list of nodes containing the AMR gene (i.e., path)
		orientation_list: list of orientation of nodes in the path
		start_pos:	the position of AMR in the first node from which the AMR gene has started
		end_pos: the position of AMR in the last node in which the AMR gene ended
	Return:
			the sequence of the AMR gene extracted from the path presented by node_list
	"""
	if not node_list:
		logging.error("ERROR: There is no node in node_list representing the AMR gene!")
		#print("ERROR: There is no node in node_list representing the AMR gene!")
		sys.exit()

	if len(node_list) == 1:
		return sequence_on_orientation(myGraph.segment(node_list[0]).sequence, orientation_list[0])[start_pos-1:end_pos]

	#Add first node
	found_amr_seq = sequence_on_orientation(myGraph.segment(node_list[0]).sequence, orientation_list[0])[start_pos-1:]
	#Add middle nodes
	if len(node_list) > 2:
		for i in range(len(node_list)-2):
			seq = sequence_on_orientation(myGraph.segment(node_list[i+1]).sequence, orientation_list[i+1])
			overlap_size = find_overlap(myGraph.segment(node_list[i]), myGraph.segment(node_list[i+1]),
										orientation_list[i], orientation_list[i+1])
			seq = seq[overlap_size:]
			found_amr_seq += seq
	#Add last node
	seq = sequence_on_orientation(myGraph.segment(node_list[-1]).sequence, orientation_list[-1])[:end_pos]
	overlap_size = find_overlap(myGraph.segment(node_list[-2]), myGraph.segment(node_list[-1]),
								orientation_list[-2], orientation_list[-1])
	seq = seq[overlap_size:]
	found_amr_seq += seq

	return found_amr_seq

# def generate_sequence_path(myGraph, node_list, orientation_list, start_pos, end_pos,
# 							pre_sequence_list, post_sequence_list,
# 							pre_path_list, post_path_list, path_dir = 'same'):
# 	"""
# 	To concatenate  the sequences and paths before and after the AMR sequence.
# 	e.g., pre_path_list = [['1-','2+']], post_path_list=[['7+'],['8+','9+','10-']],
# 	node_dir_list = ['[13+','14-]']
# 	- if path_dir = 'same' --> 		[
# 						['1-','2+','[13+','14-]','7+'],
# 						['1-','2+','[13+','14-]','8+','9+','10-']
# 									]
# 	- if path_dir = 'reverse' --> 	[
# 						['7-', '[14+', '13-]', '2-', '1+'],
# 						['10+', '9-', '8-', '[14+', '13-]', '2-', '1+']
# 									]
# 	Parameters:
# 		myGraph: The original graph
# 		node_list:	the list of nodes containing the AMR gene (i.e., path)
# 		orientation_list: list of orientation of nodes in the path
# 		start_pos:	the position of AMR in the first node from which the AMR gene has started
# 		end_pos: the position of AMR in the last node in which the AMR gene ended
# 		pre_sequence_list: the list of sequences preceding the AMR sequence in the graph
# 		post_sequence_list: the list of sequences following the AMR sequence in the graph
# 		pre_path_list: the list of nodes (and their orientation) preceding the AMR node
# 		post_path_list: the list of nodes (and their orientation) following the AMR node
# 		path_dir: 'same' if the AMR path direction (node_dir_list) is the same as what
# 			blast returned (presenting actual AMR sequence) and 'reverse' if the AMR path
# 			direction (node_dir_list) is the reverse of blast result (presenting the
# 			reverse complement of AMR sequence); in the latter, we need to return the
# 			reverse complement of the concatenated sequence to have the actual AMR
# 			sequence and not its reverse.
#
# 	"""
# 	#extract found AMR sequence
# 	found_amr_seq = extract_found_amr(myGraph, node_list, orientation_list, start_pos, end_pos)
# 	#to be able to go over both for loops make lists non-empty
# 	if not pre_path_list:
# 		pre_path_list = [[]]
# 	if not post_path_list:
# 		post_path_list = [[]]
# 	if not pre_sequence_list:
# 		pre_sequence_list = [['']]
# 	if not post_sequence_list:
# 		post_sequence_list = [['']]
# 	#add direction to nodes representing AMR gene
# 	node_orient_list = [node + orient for node, orient in zip(node_list, orientation_list)]
# 	#annotate the path representing AMR gene
# 	node_orient_list[0]='['+node_orient_list[0]
# 	node_orient_list[-1]=node_orient_list[-1]+']'
# 	#Generating all found sequences	and paths
# 	path_list = []
# 	sequence_list = []
# 	for pre_seq, pre_path in zip(pre_sequence_list, pre_path_list):
# 		for post_seq, post_path in zip(post_sequence_list, post_path_list):
# 			if path_dir == 'same':
# 				sequence_list.append(pre_seq.upper()+ found_amr_seq.lower() +post_seq.upper())
# 				path_list.append(pre_path +node_orient_list + post_path)
# 			elif path_dir == 'reverse':
# 				ssequence_list.append(rc(pre_seq.upper()+ found_amr_seq.lower() +post_seq.upper()))
# 				path_list.append(reverse_path(pre_path +node_orient_list + post_path))
#
# 	return sequence_list, path_list

def generate_sequence_path(myGraph, node_list, orientation_list, start_pos, end_pos,
							pre_sequence_list, post_sequence_list,
							pre_path_list, post_path_list, path_length_list,
							output_name, path_dir = 'same'):
	"""
	To concatenate  the sequences and paths before and after the AMR sequence.
	also, for any new sequence we add it if no similar enough sequence exists in the
	list. If such a sequence already exits in the list and its length is shorter,
	we replace it with the new sequence.
	e.g., pre_path_list = [['1-','2+']], post_path_list=[['7+'],['8+','9+','10-']],
	node_dir_list = ['[13+','14-]']
	- if path_dir = 'same' --> 		[
						['1-','2+','[13+','14-]','7+'],
						['1-','2+','[13+','14-]','8+','9+','10-']
									]
	- if path_dir = 'reverse' --> 	[
						['7-', '[14+', '13-]', '2-', '1+'],
						['10+', '9-', '8-', '[14+', '13-]', '2-', '1+']
									]
	Parameters:
		myGraph: The original graph
		node_list:	the list of nodes containing the AMR gene (i.e., path)
		orientation_list: list of orientation of nodes in the path
		start_pos:	the position of AMR in the first node from which the AMR gene has started
		end_pos: the position of AMR in the last node in which the AMR gene ended
		pre_sequence_list: the list of sequences preceding the AMR sequence in the graph
		post_sequence_list: the list of sequences following the AMR sequence in the graph
		pre_path_list: the list of nodes (and their orientation) preceding the AMR node
		post_path_list: the list of nodes (and their orientation) following the AMR node
		path_dir: 'same' if the AMR path direction (node_dir_list) is the same as what
			blast returned (presenting actual AMR sequence) and 'reverse' if the AMR path
			direction (node_dir_list) is the reverse of blast result (presenting the
			reverse complement of AMR sequence); in the latter, we need to return the
			reverse complement of the concatenated sequence to have the actual AMR
			sequence and not its reverse.

	"""
	#extract found AMR sequence
	found_amr_seq = extract_found_amr(myGraph, node_list, orientation_list, start_pos, end_pos)
	#create a temporaty directory
	temp_dir = 'temp_comparison2_'+output_name+'/'
	if not os.path.exists(temp_dir):
		try:
			os.makedirs(temp_dir)
		except OSError as exc:
			if exc.errno != errno.EEXIST:
				raise
			pass
	#to be able to go over both 'for' loops make lists non-empty
	if not pre_path_list:
		pre_path_list = [[]]
	if not post_path_list:
		post_path_list = [[]]
	if not pre_sequence_list:
		pre_sequence_list = [['']]
	if not post_sequence_list:
		post_sequence_list = [['']]
	#add direction to nodes representing AMR gene
	node_orient_list = [node + orient for node, orient in zip(node_list, orientation_list)]
	#annotate the path representing AMR gene
	node_orient_list[0]='['+node_orient_list[0]
	node_orient_list[-1]=node_orient_list[-1]+']'
	#Generating all found sequences	and paths
	path_list = []
	sequence_list = []
	counter = 0
	path_info_list = []
	for pre_seq, pre_path in zip(pre_sequence_list, pre_path_list):
		for post_seq, post_path in zip(post_sequence_list, post_path_list):
			if path_dir == 'same':
				sequence = pre_seq.upper()+ found_amr_seq.lower() +post_seq.upper()
				path = pre_path +node_orient_list + post_path
			elif path_dir == 'reverse':
				sequence = rc(pre_seq.upper()+ found_amr_seq.lower() +post_seq.upper())
				path = reverse_path(pre_path +node_orient_list + post_path)
			index, found = similar_sequence_exits(sequence_list, sequence, temp_dir)
			if not found:
				sequence_list.append(sequence)
				path_list.append(path)
				path_info_list.append(path_length_list[counter])
			elif index>=0:
				sequence_list[index] = sequence
				path_list[index] = path
				path_info_list[index] = path_length_list[counter]
			counter+=1

	# delete the temp folder
	if os.path.exists(temp_dir):
		try:
			shutil.rmtree(temp_dir)
		except OSError as e:
			logging.error("Error: %s - %s." % (e.filename, e.strerror))
			#print("Error: %s - %s." % (e.filename, e.strerror))

	return sequence_list, path_list, path_info_list

# def generate_sequence_path(myGraph, node_list, orientation_list, start_pos, end_pos,
# 							pre_sequence_list, post_sequence_list,
# 							pre_path_list, post_path_list, path_dir = 'same'):
# 	"""
# 	To concatenate  the sequences and paths before and after the AMR sequence
# 	Parameters:
# 		myGraph: The original graph
# 		node_list:	the list of nodes containing the AMR gene (i.e., path)
# 		orientation_list: list of orientation of nodes in the path
# 		start_pos:	the position of AMR in the first node from which the AMR gene has started
# 		end_pos: the position of AMR in the last node in which the AMR gene ended
# 		pre_sequence_list: the list of sequences preceding the AMR sequence in the graph
# 		post_sequence_list: the list of sequences following the AMR sequence in the graph
# 		pre_path_list: the list of nodes (and their orientation) preceding the AMR node
# 		post_path_list: the list of nodes (and their orientation) following the AMR node
# 		path_dir: 'same' if the AMR path direction (node_dir_list) is the same as what
# 			blast returned (presenting actual AMR sequence) and 'reverse' if the AMR path
# 			direction (node_dir_list) is the reverse of blast result (presenting the
# 			reverse complement of AMR sequence); in the latter, we need to return the
# 			reverse complement of the concatenated sequence to have the actual AMR
# 			sequence and not its reverse.
#
# 	"""
# 	#extract found AMR sequence
# 	found_amr_seq = extract_found_amr(myGraph, node_list, orientation_list, start_pos, end_pos)
#
# 	#create a temporaty directory
# 	temp_dir = output_dir+'temp_comparison/'
# 	if not os.path.exists(temp_dir):
# 		os.makedirs(temp_dir)
# 	#Generating all found sequences
# 	sequence_list = []
# 	to_be_removed = []
# 	for pre_seq in pre_sequence_list:
# 		for post_seq in post_sequence_list:
# 			if path_dir == 'same':
# 				sequence = pre_seq.upper()+ found_amr_seq.lower() +post_seq.upper()
# 			elif path_dir == 'reverse':
# 				sequence = rc(pre_seq.upper()+ found_amr_seq.lower() +post_seq.upper())
#
# 			sequence_list.append(sequence)
#
# 	# delete the temp folder
# 	if os.path.exists(temp_dir):
# 		try:
# 			shutil.rmtree(temp_dir)
# 		except OSError as e:
# 			print("Error: %s - %s." % (e.filename, e.strerror))
#
# 	#add direction to nodes representing AMR gene
# 	node_orient_list = [node + orient for node, orient in zip(node_list, orientation_list)]
# 	#annotate the path representing AMR gene
# 	node_orient_list[0]='['+node_orient_list[0]
# 	node_orient_list[-1]=node_orient_list[-1]+']'
# 	#Generating all found paths
# 	path_list = generate_all_paths(pre_path_list, post_path_list, node_orient_list, path_dir)
#
# 	return sequence_list, path_list

def append_path_sequence(sequence, path, sequence_list, path_list, output_dir,
							path_length, path_length_list):
	"""
	"""
	index, found = similar_sequence_exits(sequence_list, sequence, output_dir)
	if not found:
		sequence_list.append(sequence)
		path_list.append(path)
		path_length_list.append(path_length)
	elif index>=0:
		sequence_list[index] = sequence
		path_list[index] = path
		path_length_list[index] = path_length
	return sequence_list, path_list, path_length_list

def find_longest_post_neighbor(node, node_orient, node_list, both_dir):
	"""
	"""
	longest_edge = None
	max_length = 0
	for edge in (node.dovetails_L + node.dovetails_R):
		if len(edge.to_segment.sequence)>max_length and\
				edge.from_segment.name == node.name and edge.to_segment.name != node.name and\
				edge.to_segment.name not in node_list and edge.from_orient == node_orient:
			longest_edge = edge
			max_length = len(edge.to_segment.sequence)
		elif both_dir and len(edge.from_segment.sequence)>max_length and\
				edge.to_segment.name == node.name and edge.from_segment.name != node.name and\
				edge.from_segment.name not in node_list and edge.to_orient == reverse_sign(node_orient):
			longest_edge = edge
			max_length = len(edge.from_segment.sequence)
	return longest_edge

def find_longest_pre_neighbor(node, node_orient, node_list, both_dir):
	"""
	"""
	longest_edge = None
	max_length = 0
	for edge in (node.dovetails_L + node.dovetails_R):
		if len(edge.from_segment.sequence)>max_length and\
				edge.to_segment.name == node.name and edge.from_segment.name != node.name and\
				edge.from_segment.name not in node_list and edge.to_orient == node_orient:
			longest_edge = edge
			max_length = len(edge.from_segment.sequence)
		elif both_dir and len(edge.to_segment.sequence)>max_length and\
				edge.from_segment.name == node.name and edge.to_segment.name != node.name and\
				edge.to_segment.name not in node_list and edge.from_orient == reverse_sign(node_orient):
			longest_edge = edge
			max_length = len(edge.to_segment.sequence)
	return longest_edge

def extract_post_sequence_recursively_both_dir(node, node_orient, current_seq, current_path,
												length, sequence_list, path_list, node_list,
												compare_dir,length_thr, path_thr,
												current_path_length, path_length_list):
	"""
	To extract recursively the sequences following the AMR gene sequence and their paths
	In this function, we do this:
	- Find node's immediate neighbors like B
	- if there is an edge node --> B and the orient of node in that is the same as node_orient,
		B is eligible to be considered in the post_path for AMR sequence.
	- if there is an edge B --> node and the orient of node in that is reverse of
		node_orient, B with thereverse orient is eligible to be considered in the
		post_path for AMR sequence.
	Parameters:
		node:		the staring node in next paths
		node_orient:the orientation of the node in the edge this function instance was called from
		current_seq:the seq found so far
		current_path: the path found so far
		length: the remained length of the sequences around the AMR gene to be extracted
		sequence_list: the list of sequences following the AMR gene sequence
		path_list:	the list of nodes (and their orientation) following the nodes presenting the AMR gene
		node_list:	the list of nodes presenting the AMR gene
		length_thr: the threshold used for recursive pre_seq and
			post_seq until this percentage of the required length remained
			after which we just extract from the longest neighbor.
		path_thr: the threshold used for recursive pre_path and post_path
			search as long as the length of the path is less that this threshold
	Return:
		modified sequence_list and path_list
	"""
	seq = ""
	found_any_edge = False
	if length > 0 and len(current_path)<path_thr:
		for edge in (node.dovetails_L + node.dovetails_R):
			#if remained length is shorter than the threshold
			if ONLY_LONGEST_EDGE and length<=length_thr:
				import pdb; pdb.set_trace()
				longest_edge = find_longest_post_neighbor(node, node_orient, node_list, True)
				if longest_edge:
					edge = longest_edge
			#the second part of condition ensures that we are finding only nodes that are not presenting the AMR gene itself
			#Logically that shouldn't be the case unless there is a loop in the network:
			# path A -> B -> C represents the AMR; however we have an edge C -> A
			#so when we are looking for the nodes following AMR , A is selected too!!!!!
			to_segment = None
			to_orient = ''
			if edge.from_segment.name == node.name and edge.to_segment.name != node.name and\
					edge.to_segment.name not in node_list and edge.from_orient == node_orient:
				to_segment = edge.to_segment
				to_orient = edge.to_orient
			elif edge.to_segment.name == node.name and edge.from_segment.name != node.name and\
					edge.from_segment.name not in node_list and edge.to_orient == reverse_sign(node_orient):
				to_segment = edge.from_segment
				to_orient = reverse_sign(edge.from_orient)
			if to_segment and to_orient!='':
				#taking care of loops in the path
				index = exist_in_path(current_path, to_segment.name+to_orient)
				if index == -1:
					found_any_edge = True
					new_seq = sequence_on_orientation(to_segment.sequence, to_orient)
					#Remove the overlap between nodes' sequences
					overlap_size = find_overlap(edge.from_segment, edge.to_segment, edge.from_orient, edge.to_orient)
					new_seq = new_seq[overlap_size:]
					seq = current_seq + new_seq
					path=list.copy(current_path)
					path.append(str(to_segment.name)+str(to_orient))
					path_length = list.copy(current_path_length)
					path_length.append(min(len(new_seq), length))
					if len(new_seq) >= length:
						sequence_list, path_list, path_length_list =\
									append_path_sequence(seq[:len(current_seq)+length],
									path, sequence_list, path_list, compare_dir,
									path_length, path_length_list)
					else:
						sequence_list, path_list, path_length_list =\
									extract_post_sequence_recursively_both_dir(
									to_segment, to_orient, seq, path, length - len(new_seq),
									sequence_list, path_list, node_list, compare_dir,
									length_thr, path_thr,
									path_length, path_length_list)
				# graph loops are specified in { }
				elif index > -1:
					current_path[index] = '{'+current_path[index]
					current_path[-1]+='}'
			if ONLY_LONGEST_EDGE and length<=length_thr:
				#we only traverse the longest edge (already did) and not all edges
				break
	if not found_any_edge:
		sequence_list, path_list, path_length_list = append_path_sequence(current_seq,
							current_path, sequence_list, path_list, compare_dir,
							current_path_length, path_length_list)
	return sequence_list, path_list, path_length_list

def extract_post_sequence_recursively(node, node_orient, current_seq, current_path,
										length, sequence_list, path_list, node_list,
										compare_dir, length_thr, path_thr,
										current_path_length, path_length_list):
	"""
	To extract recursively the sequences following the AMR gene sequence and their paths
	Parameters:
		node:		the staring node in next paths
		node_orient:the orientation of the node in the edge this function instance was called from
		current_seq:the seq found so far
		current_path: the path found so far
		length: the remained length of the sequences around the AMR gene to be extracted
		sequence_list: the list of sequences following the AMR gene sequence
		path_list:	the list of nodes (and their orientation) following the nodes presenting the AMR gene
		node_list:	the list of nodes presenting the AMR gene
		length_thr: the threshold used for recursive pre_seq and
			post_seq until this percentage of the required length remained
			after which we just extract from the longest neighbor.
		path_thr: the threshold used for recursive pre_path and post_path
			search as long as the length of the path is less that this threshold
	Return:
		modified sequence_list and path_list

	"""
	seq = ""
	found_any_edge = False
	if length > 0 and len(current_path)<path_thr:
		for edge in (node.dovetails_L + node.dovetails_R):
			#if remained length is shorter than the threshold
			if ONLY_LONGEST_EDGE and length<=length_thr:
				longest_edge = find_longest_post_neighbor(node, node_orient, node_list, False)
				if longest_edge:
					edge = longest_edge
			#the second part of condition ensures that we are finding only nodes that are not presenting the AMR gene itself
			#Logically that shouldn't be the case unless there is a loop in the network:
			# path A -> B -> C represents the AMR; however we have an edge C -> A
			#so when we are looking for the nodes following AMR , A is selected too!!!!!
			if edge.from_segment.name == node.name and edge.to_segment.name != node.name and\
					edge.to_segment.name not in node_list and edge.from_orient == node_orient:
				#taking care of loops in the path
				index = exist_in_path(current_path, edge.to_segment.name+edge.to_orient)
				if index == -1:
					found_any_edge = True
					new_seq = sequence_on_orientation(edge.to_segment.sequence, edge.to_orient)
					#Remove the overlap between nodes' sequences
					overlap_size = find_overlap(node, edge.to_segment, node_orient, edge.to_orient)
					new_seq = new_seq[overlap_size:]
					seq = current_seq + new_seq
					path=list.copy(current_path)
					path.append(str(edge.to_segment.name)+str(edge.to_orient))
					path_length = list.copy(current_path_length)
					path_length.append(min(len(new_seq), length))
					if len(new_seq) >= length:
						sequence_list, path_list, path_length_list =\
								append_path_sequence(seq[:len(current_seq)+length],
								path, sequence_list, path_list, compare_dir,
								path_length, path_length_list)
					else:
						sequence_list, path_list, path_length_list =\
								extract_post_sequence_recursively(
								edge.to_segment, edge.to_orient, seq, path,
								length - len(new_seq), sequence_list, path_list,
								node_list, compare_dir, length_thr, path_thr,
								path_length, path_length_list)
				# graph loops are specified in { }
				elif index > -1:
					current_path[index] = '{'+current_path[index]
					current_path[-1]+='}'
			if ONLY_LONGEST_EDGE and length<=length_thr:
				#we only traverse the longest edge (already did) and not all edges
				break
	if not found_any_edge:
		sequence_list, path_list, path_length_list = append_path_sequence(current_seq, current_path,
										sequence_list, path_list, compare_dir,
										current_path_length, path_length_list)
	return sequence_list, path_list, path_length_list

def extract_post_sequence(node, node_orient, node_list, length, end_pos, compare_dir,
							remained_len_thr, path_thr, both_dir = BOTH_DIR_RECURSIVE):
	"""
	The initial function to extract the post_sequence of AMR.
	it adds the sequence following the AMR sequence (if there is any) in the last
	node in AMR path and if more is required, calls a recursive function to go over neighbors.
	Parameters:
		node:		the staring node in next paths
		node_orient:the orientation of the node in AMR path
		node_list:	the list of nodes presenting the AMR gene
		length: 	the length of the sequences around the AMR gene to be extracted
		end_pos: 	the position of AMR in the last node in which the AMR gene ended
		remained_len_thr: the threshold used for recursive pre_seq and
			post_seq until this percentage of the required length remained
			after which we just extract from the longest neighbor.
		path_thr: the threshold used for recursive pre_path and post_path
			search as long as the length of the path is less that this threshold
		both_dir:	if True we check node's orientation as well as their direction to find
		 	valid nodes in post_path/post_sequence
	Return:
		the lists of sequences and paths following AMR
	"""
	post_sequence_list = []
	post_path_list = []
	path_length = []
	path_length_list = []
	post_sequence = ''
	postfix_length = length
	if end_pos > 0:
		#attach the end of the AMR node (not included in AMR seq) to post_sequence
		postfix_length -= (len(node.sequence) - end_pos)
		post_sequence = sequence_on_orientation(node.sequence, node_orient)[end_pos:]
	#Find the sequence after the AMR gene
	if postfix_length <= 0:
		post_sequence = post_sequence[:length]
		post_sequence_list.append(post_sequence)
		path_length.append(len(post_sequence))
		path_length_list.append(path_length)
	#check all edges started from last_segment
	else:
		if len(post_sequence)>0:
			path_length.append(len(post_sequence))
		if not both_dir:
			post_sequence_list, post_path_list, path_length_list =\
									extract_post_sequence_recursively(
									node = node, node_orient = node_orient,
									current_seq = post_sequence, current_path = [],
									length = postfix_length, sequence_list = [],
									path_list = [], node_list = node_list,
									comapre_dir = compare_dir,
									length_thr = remained_len_thr, path_thr = path_thr,
									current_path_length = path_length,
									path_length_list = path_length_list)
		else:
			post_sequence_list, post_path_list, path_length_list =\
									extract_post_sequence_recursively_both_dir(
									node = node, node_orient = node_orient,
									current_seq = post_sequence, current_path = [],
									length = postfix_length, sequence_list = [],
									path_list = [], node_list = node_list,
									compare_dir = compare_dir,
									length_thr = remained_len_thr, path_thr = path_thr,
									current_path_length = path_length,
									path_length_list = path_length_list)

	return post_sequence_list, post_path_list, path_length_list

def extract_pre_sequence_recursively_both_dir(node, node_orient, current_seq, current_path,
												length, sequence_list, path_list, node_list,
												compare_dir, length_thr, path_thr,
												current_path_length, path_length_list):
	"""
	To extract recursively the sequences preceding the AMR gene sequence and their paths
	In this function, we do this:
	- Find node's immediate neighbors like B
	- if there is an edge B --> node and the orient of node in that is the same as node_orient,
		B is eligible to be considered in the pre_path for AMR sequence.
	- if there is an edge node --> B and the orient of node in that is reverse of
		node_orient, B with the reverse orient is eligible to be considered in the
		pre_path for AMR sequence.
	Parameters:
		node:		the staring node in next paths
		node_orient:the orientation of the node in the edge this function instance was called from
		current_seq:the seq found so far
		current_path: the path found so far
		length: the remained length of the sequences around the AMR gene to be extracted
		sequence_list: the list of sequences preceding the AMR gene sequence
		path_list:	the list of nodes (and their orientation) preceding the nodes presenting the AMR gene
		node_list:	the list of nodes presenting the AMR gene
		length_thr: the threshold used for recursive pre_seq and
			post_seq until this percentage of the required length remained
			after which we just extract from the longest neighbor.
		path_thr: the threshold used for recursive pre_path and post_path
			search as long as the length of the path is less that this threshold
	Return:
		modified sequence_list and path_list
	"""
	seq = ""
	found_any_edge = False
	if length > 0 and len(current_path)<path_thr:
		for edge in (node.dovetails_R + node.dovetails_L):
			#if remained length is shorter than the threshold
			if ONLY_LONGEST_EDGE and length<=length_thr:
				longest_edge = find_longest_pre_neighbor(node, node_orient, node_list, True)
				if longest_edge:
					edge = longest_edge
			#the second part ofcondition ensures that we are finding only nodes that are not presenting the AMR gene itself
			#Logically that shouldn't be the case unless there is a loop in the network:
			# path A -> B -> C represents the AMR; however we have an edge C -> A
			#so when we are looking for the nodes preceding AMR , C is selected too!!!!!
			from_segment = None
			from_orient = ''
			if edge.to_segment.name == node.name and edge.from_segment.name != node.name and\
					edge.from_segment.name not in node_list and edge.to_orient == node_orient:
				from_segment = edge.from_segment
				from_orient = edge.from_orient
			elif edge.from_segment.name == node.name and edge.to_segment.name != node.name and\
					edge.to_segment.name not in node_list and edge.from_orient == reverse_sign(node_orient):
					from_segment = edge.to_segment
					from_orient = reverse_sign(edge.to_orient)
			if from_segment and from_orient!='':
				#taking care of loops in the path
				index = exist_in_path(current_path, from_segment.name+from_orient)
				if index == -1:
					found_any_edge = True
					new_seq = sequence_on_orientation(from_segment.sequence, from_orient)
					#Remove the overlap between nodes' sequences
					overlap_size = find_overlap(edge.from_segment, edge.to_segment, edge.from_orient, edge.to_orient)
					new_seq = new_seq[:len(new_seq)-overlap_size]
					seq =  new_seq + current_seq
					path=list.copy(current_path)
					path.insert(0, str(from_segment.name) + str(from_orient))
					path_length = list.copy(current_path_length)
					path_length.insert(0, min(len(new_seq), length))
					if len(new_seq) >= length:
						sequence_list, path_list, path_length_list =\
									append_path_sequence(seq[len(new_seq)-length:],
									path, sequence_list, path_list, compare_dir,
									path_length, path_length_list)
					else:
						sequence_list, path_list, path_length_list =\
									extract_pre_sequence_recursively_both_dir(
									from_segment, from_orient, seq, path,
									length - len(new_seq), sequence_list,
									path_list, node_list, compare_dir,
									length_thr, path_thr, path_length,
									path_length_list)
				# graph loops are specified in { }
				elif index>-1:
					current_path[index] = current_path[index] + '}'
					current_path[0]='{' + current_path[0]
			if ONLY_LONGEST_EDGE and length<=length_thr:
				#we only traverse the longest edge (already did) and not all edges
				break
	if not found_any_edge:
		sequence_list, path_list, path_length_list = append_path_sequence(current_seq,
										current_path,sequence_list, path_list, compare_dir,
										current_path_length, path_length_list)
	return sequence_list, path_list, path_length_list

def extract_pre_sequence_recursively(node, node_orient, current_seq, current_path,
										length, sequence_list, path_list, node_list,
										compare_dir, length_thr, path_thr,
										current_path_length, path_length_list):
	"""
	To extract recursively the sequences preceding the AMR gene sequence and their paths
	Parameters:
		node:		the staring node in next paths
		node_orient:the orientation of the node in the edge this function instance was called from
		current_seq:the seq found so far
		current_path: the path found so far
		length: the remained length of the sequences around the AMR gene to be extracted
		sequence_list: the list of sequences preceding the AMR gene sequence
		path_list:	the list od nodes (and their orientation) preceding the nodes presenting the AMR gene
		node_list:	the list of nodes presenting the AMR gene
		length_thr: the threshold used for recursive pre_seq and
			post_seq until this percentage of the required length remained
			after which we just extract from the longest neighbor.
		path_thr: the threshold used for recursive pre_path and post_path
			search as long as the length of the path is less that this threshold
	Return:
		modified sequence_list and path_list
	"""
	seq = ""
	found_any_edge = False
	if length > 0 and len(current_path)<path_thr:
		for edge in (node.dovetails_R + node.dovetails_L):
			#if remained length is shorter than the threshold
			if ONLY_LONGEST_EDGE and length<=length_thr:
				longest_edge = find_longest_pre_neighbor(node, node_orient, node_list, False)
				if longest_edge:
					edge = longest_edge
			#the second part ofcondition ensures that we are finding only nodes that are not presenting the AMR gene itself
			#Logically that shouldn't be the case unless there is a loop in the network:
			# path A -> B -> C represents the AMR; however we have an edge C -> A
			#so when we are looking for the nodes preceding AMR , C is selected too!!!!!
			if edge.to_segment.name == node.name and edge.from_segment.name != node.name and\
					edge.from_segment.name not in node_list and edge.to_orient == node_orient:
				#taking care of loops in the path
				index = exist_in_path(current_path, edge.from_segment.name+edge.from_orient)
				if index == -1:
					found_any_edge = True
					new_seq = sequence_on_orientation(edge.from_segment.sequence, edge.from_orient)
					#Remove the overlap between nodes' sequences
					overlap_size = find_overlap(edge.from_segment, node, edge.from_orient, node_orient)
					new_seq = new_seq[:len(new_seq)-overlap_size]
					seq =  new_seq + current_seq
					path=list.copy(current_path)
					path.insert(0, str(edge.from_segment.name) + str(edge.from_orient))
					path_length = list.copy(current_path_length)
					path_length.insert(0, min(len(new_seq), length))
					if len(new_seq) >= length:
						sequence_list, path_list, path_length_list =\
									append_path_sequence(seq[len(new_seq)-length:],
									path, sequence_list, path_list, compare_dir,
									path_length, path_length_list)
					else:
						sequence_list, path_list, path_length_list =\
						 			extract_pre_sequence_recursively(
									edge.from_segment, edge.from_orient, seq, path,
									length - len(new_seq), sequence_list, path_list,
									node_list, compare_dir, length_thr, path_thr,
									path_length, path_length_list)
				# graph loops are specified in { }
				elif index>-1:
					current_path[index] = current_path[index] + '}'
					current_path[0]='{' + current_path[0]
			if ONLY_LONGEST_EDGE and length<=length_thr:
				#we only traverse the longest edge (already did) and not all edges
				break
	if not found_any_edge:
		sequence_list, path_list, path_length_list = append_path_sequence(current_seq,
									current_path, sequence_list, path_list,
									compare_dir, current_path_length, path_length_list)
	return sequence_list, path_list, path_length_list

def extract_pre_sequence(node, node_orient, node_list, length, start_pos, compare_dir,
							remained_len_thr, path_thr, both_dir = BOTH_DIR_RECURSIVE):
	"""
	The initial function to extract the pre_sequence of AMR.
	it adds the sequence preceding the AMR sequence (if there is any) in the first
	node in AMR path and if more sequence is required, calls a recursive function to go over neighbors.
	Parameters:
		node:		the staring node in next paths
		node_orient:the orientation of the node in AMR path
		node_list:	the list of nodes presenting the AMR gene
		length: 	the length of the sequences around the AMR gene to be extracted
		end_pos: 	the position of AMR in the last node in which the AMR gene ended
		remained_len_thr: the threshold used for recursive pre_seq and
			post_seq until this percentage of the required length remained
			after which we just extract from the longest neighbor.
		path_thr: the threshold used for recursive pre_path and post_path
			search as long as the length of the path is less that this threshold
		both_dir:	if True we check node's orientation as well as their direction to find
		 	valid nodes in pre_path/pre_sequence
	Return:
		the lists of sequences and paths preceding AMR

	"""
	pre_sequence_list = []
	pre_path_list = []
	pre_sequence = ''
	path_length = []
	path_length_list = []
	prefix_length = length
	if start_pos > 0:
		#attach the beginning of the AMR node (not included in AMR seq) to pre_sequence
		prefix_length -= (start_pos-1)
		pre_sequence = sequence_on_orientation(node.sequence, node_orient)[:start_pos-1]
	#Find the sequence before the AMR gene
	if prefix_length <= 0:
		pre_sequence = pre_sequence[-length:]
		pre_sequence_list.append(pre_sequence)
		path_length.append(len(pre_sequence))
		path_length_list.append(path_length)
	#check all edges ended at first node
	else:
		if len(pre_sequence)>0:
			path_length.append(len(pre_sequence))
		if not both_dir:
			pre_sequence_list, pre_path_list, path_length_list =\
									extract_pre_sequence_recursively(
									node = node, node_orient = node_orient,
									current_seq = pre_sequence, current_path = [],
									length = prefix_length, sequence_list = [],
									path_list = [], node_list = node_list,
									compare_dir = compare_dir,
									length_thr = remained_len_thr,path_thr = path_thr,
									current_path_length = path_length,
									path_length_list = path_length_list)
		else:
			pre_sequence_list, pre_path_list, path_length_list =\
									extract_pre_sequence_recursively_both_dir(
									node = node, node_orient = node_orient,
									current_seq = pre_sequence, current_path = [],
									length = prefix_length, sequence_list = [],
									path_list = [], node_list = node_list,
									compare_dir = compare_dir,
									length_thr = remained_len_thr, path_thr = path_thr,
									current_path_length = path_length,
									path_length_list = path_length_list)

	return pre_sequence_list, pre_path_list, path_length_list

def find_reverse_blast_path(myGraph, node_list, orientation_list, start_pos, end_pos):
	"""
	To reverse the AMR path returned by blast and also adjust start_pos and end_pos
	e.g., [A+, B-, C-] --> [C+, B+, A-]
							new_node_list = [C, B, A]
							new_orientation_list = [+, +, -]
							new_start_pos = len(C) - end_pos + 1
							new_end_pos = len(A) -  start_pos +1
	Parameters:
		myGraph: The assembly graph
		node_list:	the list of nodes containing the AMR gene (i.e., path)
		orientation_list: list of orientation of nodes in the AMR path
		start_pos:	the position of AMR in the first node from which the AMR gene has started
		end_pos: the position of AMR in the last node in which the AMR gene ended
	Return:
	"""
	node_list_reverse = list.copy(node_list)
	orientation_list_reverse = list.copy(orientation_list)
	#single node in AMR path
	if len(node_list)==1:
		node = myGraph.segment(node_list[0])
		orientation_list_reverse[0] = reverse_sign(orientation_list[0])
		start_pos_reverse = len(node.sequence)-end_pos+1
		end_pos_reverse = len(node.sequence)-start_pos+1
	#multiple nodes in AMR path
	elif len(node_list)>1:
		node_list_reverse.reverse()
		orientation_list_reverse.reverse()
		for i in range(len(orientation_list_reverse)):
			orientation_list_reverse[i]=reverse_sign(orientation_list_reverse[i])
		node1 = myGraph.segment(node_list_reverse[0])
		node2 = myGraph.segment(node_list_reverse[-1])
		start_pos_reverse = len(node1.sequence)-end_pos+1
		end_pos_reverse = len(node2.sequence)-start_pos+1

	return node_list_reverse, orientation_list_reverse, start_pos_reverse, end_pos_reverse

def write_paths_info_to_file(paths_info_list, paths_info_file, seq_counter):
	"""
	"""
	counter = seq_counter
	with open(paths_info_file,'a') as fd:
		writer = csv.writer(fd)
		for i, path_info in enumerate(paths_info_list):
			for node_info in path_info:
				counter = i+seq_counter+1
				writer.writerow([counter, node_info['node'], node_info['coverage'],
									node_info['start'], node_info['end']])
	return counter


def generate_node_range_coverage(myGraph, node_list, orientation_list, start_pos,
								end_pos, pre_path_list, pre_path_length_list,
								post_path_list, post_path_length_list,
								max_kmer_size, assembler = Assembler_name.meta_spades):
	"""
	"""
	if not pre_path_list:
		pre_path_list = [[]]
	if not post_path_list:
		post_path_list = [[]]
	#find amr_path_length_info
	start = 0
	amr_path_length_info = []
	for i, node in enumerate(node_list):
		segment = myGraph.segment(node)
		if assembler == Assembler_name.meta_spades or assembler == Assembler_name.bcalm:
			coverage = segment.KC/(len(segment.sequence)-max_kmer_size)
		elif assembler == Assembler_name.megahit:
			coverage = float(node.split('cov_')[1].split('_')[0])
		else:
			logging.error("no way of calculating node coverage has been defined for this assembler!")
			sys.exit()
		#the length of first node in amr
		if i == 0 and len(node_list)==1:
			length = end_pos - start_pos + 1
		elif i==0 and len(node_list)>1:
			length = len(segment.sequence) -  start_pos + 1
		#the length of intermediate nodes in amr
		elif i < (len(node_list) - 1):
			overlap_size = find_overlap(myGraph.segment(node_list[i-1]), segment,
										orientation_list[i-1], orientation_list[i])
			length = len(segment.sequence) - overlap_size
		#the length of last node in amr
		else:
			overlap_size = find_overlap(myGraph.segment(node_list[i-1]), segment,
										orientation_list[i-1], orientation_list[i])
			length = end_pos - overlap_size
		node_info = {'node': node, 'coverage': coverage, 'start': start, 'end': start+length-1}
		amr_path_length_info.append(node_info)
		start =  start+length

	#pre_path_length_info
	pre_paths_info =[]
	for path, path_length in zip(pre_path_list, pre_path_length_list):
		start = 0
		pre_path_info=[]
		assert (len(path)==len(path_length) or len(path)==len(path_length) - 1),"inconsistent length of arrays: path vs path_length"
		difference = len(path_length) - len(path)
		end_common = len(path_length) -  difference
		for node, length in zip(path, path_length[:end_common]):
			pure_node = find_node_name(node)
			segment = myGraph.segment(pure_node)
			if assembler == Assembler_name.meta_spades or assembler == Assembler_name.bcalm:
				coverage = segment.KC/(len(segment.sequence)-max_kmer_size)
			elif assembler == Assembler_name.megahit:
				coverage = float(pure_node.split('cov_')[1].split('_')[0])
			else:
				logging.error("no way of calculating node coverage has been defined for this assembler!")
				sys.exit()
			node_info = {'node': pure_node, 'coverage': coverage, 'start': start, 'end': start+length-1}
			start = start+length
			pre_path_info.append(node_info)
		if len(path)==len(path_length) - 1:
			segment = myGraph.segment(node_list[0])
			if assembler == Assembler_name.meta_spades or assembler == Assembler_name.bcalm:
				coverage = segment.KC/(len(segment.sequence)-max_kmer_size)
			elif assembler == Assembler_name.megahit:
				coverage = float(node_list[0].split('cov_')[1].split('_')[0])
			else:
				logging.error("no way of calculating node coverage has been defined for this assembler!")
				sys.exit()
			node_info_last = {'node': node_list[0], 'coverage': coverage, 'start': start, 'end': start+path_length[-1]-1}
			pre_path_info.append(node_info_last)
		pre_paths_info.append(pre_path_info)

	#attach amr info to pre_path
	for i in range(len(pre_paths_info)):
		if len(pre_paths_info[i])>0:
			lag = pre_paths_info[i][-1]['end']+1
		else:
			lag=0
		for amr_info in amr_path_length_info:
			tmp = amr_info.copy()
			tmp['start'] += lag
			tmp['end'] += lag
			pre_paths_info[i].append(tmp)

	#post_path_length_info
	post_paths_info =[]
	for path, path_length in zip(post_path_list, post_path_length_list):
		start = 0
		post_path_info=[]
		assert (len(path)==len(path_length) or len(path)==len(path_length) - 1),"inconsistent length of arrays: path vs path_length"
		start_index = len(path_length) - len(path)
		if len(path)==len(path_length) - 1:
			segment = myGraph.segment(node_list[-1])
			if assembler == Assembler_name.meta_spades or assembler == Assembler_name.bcalm:
				coverage = segment.KC/(len(segment.sequence)-max_kmer_size)
			elif assembler == Assembler_name.megahit:
				coverage = float(node_list[-1].split('cov_')[1].split('_')[0])
			else:
				logging.error("no way of calculating node coverage has been defined for this assembler!")
				sys.exit()
			node_info = {'node': node_list[-1], 'coverage': coverage, 'start': start, 'end': start+path_length[0]-1}
			post_path_info.append(node_info)
			start = start + path_length[0]
		for node, length in zip(path, path_length[start_index:]):
			pure_node = find_node_name(node)
			segment = myGraph.segment(pure_node)
			if assembler == Assembler_name.meta_spades or assembler == Assembler_name.bcalm:
				coverage = segment.KC/(len(segment.sequence)-max_kmer_size)
			elif assembler == Assembler_name.megahit:
				coverage = float(pure_node.split('cov_')[1].split('_')[0])
			else:
				logging.error("no way of calculating node coverage has been defined for this assembler!")
				sys.exit()
			node_info = {'node': pure_node, 'coverage': coverage, 'start': start, 'end': start+length-1}
			start = start+length
			post_path_info.append(node_info)
		post_paths_info.append(post_path_info)

	#generate the entire path info by adding post_path
	paths_info = []
	for pre_path in pre_paths_info:
		if len(pre_path)>0:
			lag = pre_path[-1]['end']+1
		else:
			lag=0
		for post_path in post_paths_info:
			path_info = list.copy(pre_path)
			for node_info in post_path:
				tmp = node_info.copy()
				tmp['start']+=lag
				tmp['end']+=lag
				path_info.append(tmp)
			paths_info.append(path_info)

	return paths_info

def extract_neighborhood_sequence(gfa_file, length, node_list, orientation_list,
									start_pos, end_pos, remained_len_thr, path_thr,
									output_name, max_kmer_size,
									assembler = Assembler_name.meta_spades):
	"""
	To extract all linear sequences with length = length from the start and end of the AMR gene
	Parameters:
		myGraph: The assembly graph
		length: the length of all sequences around the AMR gene to be extracted
		node_list:	the list of nodes containing the AMR gene (i.e., path)
		orientation_list: list of orientation of nodes in the path
		start_pos:	the position of AMR in the first node from which the AMR gene has started
		end_pos: the position of AMR in the last node in which the AMR gene ended
		remained_len_thr: the threshold used for recursive pre_seq and
			post_seq until this percentage of the required length remained
			after which we just extract from the longest neighbor.
		path_thr: the threshold used for recursive pre_path and post_path
			search as long as the length of the path is less that this threshold
	"""
	myGraph = gfapy.Gfa.from_file(gfa_file)

	last_segment = myGraph.segment(node_list[-1])
	if start_pos == 0:
		start_pos = 1
	if end_pos == 0:
		end_pos = len(str(last_segment.sequence))

	logging.debug('last_segment = '+last_segment.name+' start_pos = '+str(start_pos)+' end_pos= '+str(end_pos))
	#create a temporaty directory
	compare_dir = 'temp_comparison_'+output_name+'/'
	if not os.path.exists(compare_dir):
		try:
			os.makedirs(compare_dir)
		except OSError as exc:
			if exc.errno != errno.EEXIST:
				raise
			pass

	#Find the sequences and paths for the path provided by bandage+blast for AMR gene
	#Find the sequence before the AMR gene
	logging.debug('Running extract_pre_sequence '+output_name+' ...')
	pre_sequence_list, pre_path_list, pre_path_length_list =\
	 							extract_pre_sequence(myGraph.segment(node_list[0]),
								orientation_list[0], node_list, length, start_pos,
								compare_dir, remained_len_thr, path_thr)
	#Find the sequence after the AMR gene
	logging.debug('Running extract_post_sequence '+output_name+' ...')
	post_sequence_list, post_path_list, post_path_length_list =\
								extract_post_sequence(last_segment,
								orientation_list[-1], node_list, length, end_pos,
								compare_dir, remained_len_thr, path_thr)
	#combine path_info from pre and post
	logging.debug('Running generate_node_range_coverage '+output_name+' ...')
	path_length_list = generate_node_range_coverage(myGraph, node_list, orientation_list, start_pos,
							end_pos, pre_path_list, pre_path_length_list,
							post_path_list, post_path_length_list, max_kmer_size, assembler)
	#Combine pre_ and post_ sequences and paths
	logging.debug('Running generate_sequence_path '+output_name+' ...')
	sequence_list, path_list, path_info_list = generate_sequence_path(myGraph,
										node_list, orientation_list, start_pos,
										end_pos, pre_sequence_list,
										post_sequence_list, pre_path_list,
										post_path_list, path_length_list,
										output_name, 'same')
	# sequence_list, path_list = generate_sequence_path(myGraph,node_list,
	# 									orientation_list, start_pos, end_pos,
	# 									pre_sequence_list, post_sequence_list,
	# 									pre_path_list, post_path_list, 'same')

	if not BOTH_DIR_RECURSIVE:
		#Find the sequences and paths for the reverse of the path provided by bandage+blast for AMR gene
		node_list_reverse, orientation_list_reverse, start_pos_reverse, end_pos_reverse =\
			find_reverse_blast_path(myGraph, node_list, orientation_list, start_pos, end_pos)
		#Find the sequence before the AMR gene
		pre_sequence_list, pre_path_list, pre_path_length_list =\
								extract_pre_sequence(myGraph.segment(node_list_reverse[0]),
								orientation_list_reverse[0], node_list_reverse, length, start_pos_reverse,
								compare_dir, remained_len_thr, path_thr)
		#Find the sequence after the AMR gene
		post_sequence_list, post_path_list, post_path_length_list =\
								extract_post_sequence(myGraph.segment(node_list_reverse[-1]),
								orientation_list_reverse[-1], node_list_reverse, length, end_pos_reverse,
								compare_dir, remained_len_thr, path_thr)
		#Combine pre_ and post_ sequences and paths
		sequence_list_reverse, path_list_reverse = generate_sequence_path(myGraph,node_list_reverse,
											orientation_list_reverse, start_pos_reverse,
											end_pos_reverse, pre_sequence_list,
											post_sequence_list, pre_path_list, post_path_list,
											output_name, 'reverse')

		sequence_list, path_list = remove_identical_sub_sequences(sequence_list+sequence_list_reverse,
																	path_list+path_list_reverse)
	# delete the temp folder
	if os.path.exists(compare_dir):
		try:
			shutil.rmtree(compare_dir)
		except OSError as e:
			logging.error("Error: %s - %s." % (e.filename, e.strerror))
			#print("Error: %s - %s." % (e.filename, e.strerror))

	return sequence_list, path_list, path_info_list

def find_amr_related_nodes(amr_file, gfa_file, output_dir, bandage_path = BANDAGE_PATH,
							threshold =  95, output_pre = '', align_file = ''):
	"""
	Run bandage+blast to find the sequence in amr_file (as the query) in the assembly
	graph (gfa_file), and extract the path(s) representing the AMR sequence (i.e.,
	AMR path) in the assembly graph.
	e.g., (121) A+, B-, C+ (143) --> paths_info =	[
										{nodes:['A', 'B', 'C'],
										orientations: ['+', '-', '+'],
										start_pos: 121,
										end_pos: 143}
													]
	Parameters:
		amr_file: the FASTA file containing the AMR sequence (our query)
		gfa_file: the GFA file containing the assembly graph (our DB)
		output_dir: the output directory to store files
		bandage_path: the address of bandage executation file
		threshold: threshold for identity and coverage
		output_pre: used for naming output file
	Return:
		A boolean value which is True if any path was found and
		A list of dictionaries each denoting an AMR path
	"""
	if align_file =='':
		#Run bandage+blast
		output_name=output_dir+output_pre+'_align_'+datetime.datetime.now().strftime('%Y-%m-%d_%H-%M-%S')
		if os.path.isfile(output_name+'.tsv'):
			os.remove(output_name+'.tsv')
		command = bandage_path +' querypaths '+gfa_file+' '+amr_file+' '+output_name
		os.system(command)
		align_file = output_name+".tsv"
	#Process the output tsv file
	paths_info = []
	found = False
	logging.debug("Reading align_file = "+align_file+" ...")
	with open(align_file) as tsvfile:
		reader = csv.reader(tsvfile, delimiter='\t')
		#skip the header
		next(reader)
		for row in reader:
			coverage = float(re.sub('[%]','',row[3]))
			identity = float(re.sub('[%]','',row[5]))
			if int(coverage) >= threshold and int(identity)>=threshold:
				#print('AMR '+amr_name+' was found!')
				found = True
				cell_info = row[1].strip()
				nodes, orientation_list, start_pos, end_pos = extract_nodes_in_path(cell_info)
				path_info = {'nodes':nodes, 'orientations':orientation_list,
								'start_pos':start_pos, 'end_pos':end_pos}
				paths_info.append(path_info)
	# if not found:
	# 	os.remove(output_name+'.tsv')
	return found, paths_info

def neighborhood_sequence_extraction(gfa_file, length, output_dir,
									bandage = BANDAGE_PATH, threshold = 95,
									seq_name_prefix = 'ng_sequences_',
									#output_name = 'ng_sequences',
									path_node_threshold = 10 ,
									path_seq_len_percent_threshod = 90,
									max_kmer_size = 55,
									assembler = Assembler_name.meta_spades,
									amr_seq_align_file = ''):
	"""
	The core function to extract the sequences/paths preceding and following the AMR sequence
	in the assembly graph
	Parameters:
		amr_file:	the FASTA file containing the AMR sequence
		gfa_file:	the GFA file containing the assembly graph
		length:		the length of all sequences around the AMR gene to be extracted
		output_dir: the output directory to store files
		bandage: 	the address of bandage executation file
		threshold: 	threshold for identity and coverage
		output_name: used for naming output file
		path_node_threshold: the threshold used for recursive pre_path and post_path
			search as long as the length of the path is less that this threshold
		path_seq_len_percent_threshod: the threshold used for recursive pre_seq and
			post_seq until we have this percentage of the required length
			after which we just extract from the longest neighbor.
	Return:
		the name of file containing the list of extracted sequences/paths
	"""
	amr_file, align_file = amr_seq_align_file
	logging.debug('amr_file = '+amr_file+'\t align_file = '+align_file)
	output_name = seq_name_prefix +os.path.splitext(os.path.basename(amr_file))[0]
	if not os.path.exists(output_dir+'alignment_files/'):
		try:
			os.makedirs(output_dir+'alignment_files/')
		except OSError as exc:
			if exc.errno != errno.EEXIST:
				raise
			pass
	remained_len_thr = length - (length*path_seq_len_percent_threshod/100.0)
	#find all AMR paths in the assembly graph
	_, amr_paths_info = find_amr_related_nodes(amr_file, gfa_file, output_dir+'alignment_files/',
											bandage, threshold, output_name, align_file)
	seq_output_dir = output_dir+'sequences/'
	if not os.path.exists(seq_output_dir):
		try:
			os.makedirs(seq_output_dir)
		except OSError as exc:
			if exc.errno != errno.EEXIST:
				raise
			pass
	seq_file = seq_output_dir+output_name+'_'+str(length)+'_'+\
	datetime.datetime.now().strftime('%Y-%m-%d_%H-%M')+'.txt'
	#csv file for path_info
	if not os.path.exists(output_dir+'paths_info/'):
		try:
			os.makedirs(output_dir+'paths_info/')
		except OSError as exc:
			if exc.errno != errno.EEXIST:
				raise
			pass
	paths_info_file = output_dir+'paths_info/'+output_name+'_'+str(length)+'_'+\
	datetime.datetime.now().strftime('%Y-%m-%d_%H-%M')+'.csv'
	with open(paths_info_file,'a') as fd:
		writer = csv.writer(fd)
		writer.writerow(['sequence', 'node', 'coverage', 'start', 'end'])
	logging.debug("Calling extract_neighborhood_sequence for "+os.path.basename(amr_file)+" ...")
	#Extract the sequenc of AMR neighborhood
	seq_counter = 0
	for i, amr_path_info in enumerate(amr_paths_info):
		node_list = amr_path_info['nodes']
		orientation_list = amr_path_info['orientations']
		start_pos = amr_path_info['start_pos']
		end_pos = amr_path_info['end_pos']
		sequence_list, path_list, path_info_list =\
							extract_neighborhood_sequence(gfa_file, length,
							node_list, orientation_list, start_pos, end_pos,
							remained_len_thr, path_node_threshold, output_name,
							max_kmer_size, assembler)
		if not sequence_list:
			continue
		write_sequences_to_file(sequence_list, path_list, seq_file)
		logging.info("NOTE: The list of neighborhood sequences (extracted from assembly graph)\
	 		has been stroed in " + seq_file)
		seq_counter = write_paths_info_to_file(path_info_list, paths_info_file, seq_counter)
	return seq_file, paths_info_file

def neighborhood_graph_extraction(gfa_file, distance, output_dir,
								bandage = BANDAGE_PATH, threshold = 95,
								seq_name_prefix = 'ng_subgraph',
								amr_seq_align_file = ''):
	"""
	"""
	amr_file, align_file = amr_seq_align_file
	logging.debug('amr_file = '+amr_file+'\t align_file = '+align_file)
	output_name = seq_name_prefix +os.path.splitext(os.path.basename(amr_file))[0]
	#Retrieve the AMR gene sequence
	amr_seq = retrieve_AMR(amr_file)
	if not os.path.exists(output_dir+'alignment_files/'):
		try:
			os.makedirs(output_dir+'alignment_files/')
		except OSError as exc:
			if exc.errno != errno.EEXIST:
				raise
			pass
	_, paths_info = find_amr_related_nodes(amr_file, gfa_file, output_dir+'alignment_files/',
											bandage, threshold, output_name, align_file)

	#Find neighbourhood of found nodes AND create a new graph contating the target subgraph
	subgraph_file_list =[]
	for i, path_info in enumerate(paths_info):
		node_list = path_info['nodes']
		#Find the list of segments and edges
		myGraph, segment_list, edge_list =extract_subgraph(node_list, distance, gfa_file)
		#Generate a graph with the lists of segments and edges
		subGraph = create_graph(myGraph, segment_list, edge_list)
		#write the subgraph in a GFA file
		subgraph_file = output_dir+output_name+'_'+str(distance)+'_'+str(i+1)+'_'+\
			datetime.datetime.now().strftime('%Y-%m-%d_%H-%M')+'.gfa'
		subGraph.to_file(subgraph_file)
		logging.info("NOTE: The extracted subgraph has been stored in " + subgraph_file)
		subgraph_file_list.append(subgraph_file)

	return subgraph_file_list

if __name__=="__main__":

	text = 'This code is used to find the context of a given AMR gene'
	parser = argparse.ArgumentParser(description=text)
	parser.add_argument('--amr','-A', type=str, default = '',
		help = 'the path of the file containing the AMR gene sequence')
	parser.add_argument('--gfa', '-G', type=str, default='',
		help = 'the GFA file containing the graph')
	parser.add_argument('--distance', '-D', type = int, default=1,
		help = 'the maximum distance of neighborhood nodes to be extracted from the AMR gene')
	parser.add_argument('--length', '-L', type = int, default=1000,
		help = 'the length of AMR gene\'s neighbourhood to be extracted')
	parser.add_argument('--output_dir','-O', type=str, default = OUT_DIR,
		help = 'the output directory to store the results')
	args = parser.parse_args()

	if not os.path.exists(args.output_dir):
		os.makedirs(args.output_dir)

	if not args.amr:
		logging.error('please enter the path for the AMR gene file')
		#print('please enter the path for the AMR gene file')
		sys.exit()

	if not args.gfa:
		logging.error('please enter the path for the assembly file')
		#print('please enter the path for the assembly file')
		sys.exit()

	neighborhood_graph_extraction(arg.amr,args.gfa, args.distance, args.output_dir)
	neighborhood_sequence_extraction(args.gfa, args.length, args.output_dir, amr_file = (arg.amr,''))
