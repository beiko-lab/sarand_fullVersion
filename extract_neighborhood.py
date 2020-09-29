"""
File:		extract_neighborhood.py
Aythor:		Somayeh Kafaie
Date:		July 2020
Purpose:	To extract the neighborhood of an AMR gene from an assembly graph

To run:
- The Bandage+BLAST implementation (type = 1):
	python extract_neighborhood.py --amr/-A <AMR gene file path in FASTA format>
	--gfa/-G <GFA assembly graph> --distance/-D <max distance from AMR gene (default=1)>
	--length/-L <lebgth of the linear sequence around AMR gene to be extracted (default = 1000)>

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
from gfapy.sequence import rc

BANDAGE_PATH = '/media/Data/tools/Bandage_Ubuntu_dynamic_v0_8_1/Bandage'
MAX_PATH_NODES = 6
OUT_DIR = 'output/'

def retrieve_AMR(file_path):
	"""
	To read the AMR gene from the text file.
	Parameters:
		file_path:	the address of the file containing the AMR gene
	Return:
		the sequence of the AMR gene in lower case
	"""
	with open(file_path) as fp:
		#for i, line in enumerate(fp):
		for line in fp:
			if line.startswith('>'):
				continue
			return line.lower()

def extract_amr_neighborhood_in_ref_genome(amr_seq, ref_path, neighborhood_len):
	"""
	To find the neighborhood of amr_seq in the genome ref_path
	Parameters:
		amr_seq:	the sequence of the AMR gene
		ref_path:	the address of the reference genome
		neighborhood_len: 	the length of neighborhood from each side to be extracted
	Return:
		the neighborhood of amr_seq with maximum length = 2*len + length(amr_seq)
			and the start position of AMR in the extracted sequence
	"""
	seq = ''
	amr_start = 0
	for record in SeqIO.parse(open(ref_path,'r'),'fasta'):
		i = str(record.seq).lower().find(amr_seq)
		if i >= 0:
			if i-neighborhood_len >= 0:
				seq = str(record.seq)[i-neighborhood_len:i+len(amr_seq)]
				amr_start = neighborhood_len
			else:
				seq = str(record.seq)[:i+len(amr_seq)]
				amr_start = i
			if (i+len(amr_seq)+neighborhood_len)<=len(str(record.seq)):
				seq+=str(record.seq)[i+len(amr_seq):i+len(amr_seq)+neighborhood_len]
			else:
				seq+=str(record.seq)[i+len(amr_seq):]
	#added by 1 because the index starts from 0 in the string
	return seq, amr_start+1

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
			print("WARNING: the link is a loop for node#"+segment.name)
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
			print("WARNING: the link is a loop for node#"+segment.name)
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
	for each node to do the job
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
	#for node_list in node_lists:
	#nodes = node_list.split(',')
	#for node in nodes:
	#for now we just find the neighborhood of the first and last nodes in the path
	#node = re.sub('[+-]', '', node)
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
		direction_list: list of direction of nodes -> e.g., [-, +]
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
	direction_list = []
	nodes = path.split(',')
	for node in nodes:
		if '-' in node:
			direction_list.append('-')
		else:
			direction_list.append('+')
		node = re.sub('[+-]', '', node.split()[0])
		node_list.append(node)
	return node_list, direction_list, start_pos, end_pos

def write_sequences_to_file(sequence_list, path_list, node_list, length, file_name):
	"""
	To write in a filethe extracted sequences and paths in neighborhood of the AMR gene
	Parameters:
		sequence_list: the list of sequences extracted in AMR gene neighborhood
		path_list: the list of all paths (a list of nodes) in AMR gene neighborhood
		node_list: the list of nodes presenting the AMR gene
		length: the lebgth of sequences extracted follwoing and preceding the AMR sequence
		file_name: the name of file to be writtn in
	"""
	file = open(file_name, 'w')
	#myLine = "> The list of paths and sequences (with length "+ str(length) +" ) around the AMR gene:\n"
	#file.write(myLine)
	for seq, path in zip(sequence_list, path_list):
		myLine = "Path-> " + ", ".join(path) + ":"
		file.write(myLine+"\n")
		file.write(seq+"\n")
	file.close()

def generate_all_paths(pre_path_list, post_path_list, node_dir_list):
	"""
	Generate the list of all paths which are combination of the input lists
	e.g., pre_path_list = [['1-','2+']], post_path_list=[['7+'],['8+','9+','10-']],
	node_dir_list = ['[13+','14-]'] --> [
						['1-','2+','[13+','14-]','7+'],
						['1-','2+','[13+','14-]','8+','9+','10-']
										]
	Parameters:
		pre_path_list: the list of nodes (and their direction) precending the AMR gene
		post_path_list: the list of nodes (and their direction) following the AMR gene
		node_dir_list: the list of nodes (and their direction) presenting the AMR gene
	"""
	path_list = []
	if pre_path_list and post_path_list:
		for pre_path in pre_path_list:
			for post_path in post_path_list:
				path_list.append(pre_path +node_dir_list + post_path)
				#path_list.append(list(pre_path +node_dir_list + post_path))
	elif pre_path_list:
		for pre_path in pre_path_list:
			path_list.append(pre_path +node_dir_list)
	elif post_path_list:
		for post_path in post_path_list:
			path_list.append(node_dir_list + post_path)
	else:
		path_list.append(node_dir_list)
	return path_list

def sequence_on_direction(seq, dir):
	"""
	To	return the sequence based on the direction used in the graph
	Parameters:
	 	dir:	+ / -
		seq:	the original sequence
	Return: reversed complement for -; otherwise seq itself
		for dir = '-' and seq = 'AACTG' -> 'CAGTT'
	"""
	if dir=='-':
			return rc(seq)
	else:
		return seq

def find_overlap(from_node, to_node, validate_overlap_size = True):
	"""
	To return the size of overlap beween  from_node and to_node
	e.g.; from_node = 2, to_node = 3, edge = [L 2 + 3 + 55M]
				---->  return 55
	Parameters:
		from_node:	the first node of the edge
		to_node:	the second node of the edge
		validate_overlap_size: to verify if the extracted size for overlap is correct and matches
	Return:
		the size of overlap
	"""
	size = 0
	for edge in (from_node.dovetails_L + from_node.dovetails_R):
		if edge.from_segment == from_node and edge.to_segment == to_node:
			if edge.overlap[0].code == 'M':
				size = edge.overlap[0].length
				if (validate_overlap_size):
					seq1 = sequence_on_direction(edge.from_segment.sequence, edge.from_orient)
					seq2 = sequence_on_direction(edge.to_segment.sequence, edge.to_orient)
					#if seq1[len(seq1)-size:] == seq2[:size]:
					#	print("The overlap matches!")
					#else:
					if seq1[len(seq1)-size:] != seq2[:size]:
						print("WARNING: the overlap doesn't match: "+str(edge))
						print("first sequence: "+seq1)
						print("second sequence: "+seq2)
			else:
				print("WARNING: the sequence from node "+from_node.name+" to node "+
						to_node.name + "has "+ edge.overlap[0] + " instead of overlap!")
			return size

	#In case for any wierd technical issue the edge wan't found using from_node--->
	for edge in (to_node.dovetails_R + to_node.dovetails_L):
		if edge.from_segment == from_node and edge.to_segment == to_node:
			if edge.overlap[0].code == 'M':
				size = edge.overlap[0].length
				if (validate_overlap_size):
					seq1 = sequence_on_direction(edge.from_segment.sequence, edge.from_orient)
					seq2 = sequence_on_direction(edge.to_segment.sequence, edge.to_orient)
					if seq1[len(seq1)-size:] == seq2[:size]:
						print("The overlap matches!")
					else:
						print("WARNING: the overlap doesn't match: "+str(edge))
						print("first sequence: "+seq1)
						print("second sequence: "+seq2)
			else:
				print("WARNING: the sequence from node "+from_node.name+" to node "+
						to_node.name + "has "+ edge.overlap[0] + " instead of overlap!")
			return size


def extract_found_amr(myGraph, node_list, direction_list, start_pos, end_pos):
	"""
	To extract the sequence representing AMR gene from the path fund by
	Bandage+BLAST presented in node_list; I could have used the original AMR gene
	sequence but I thought there might be cases that the assembler can't assemble
	the entire AMR gene or there are some errors in assembled one.
	Parameters:
		myGraph: The original graph
		node_list:	the list of nodes containing the AMR gene (i.e., path)
		direction_list: list of direction of nodes in the path
		start_pos:	the position of the first node from which the AMR gene has started
		end_pos: the position of the last node in which the AMR gene ended
	Return:
			the sequence of the AMR gene extracted from the path presented by node_list
	"""
	if start_pos == 0:
		start_pos = 1
	if end_pos == 0:
		end_pos = len(myGraph.segment(node_list[-1]).sequence)

	if not node_list:
		print("ERROR: There is no node in node_list representing the AMR gene!")
		sys.exit()

	if len(node_list) == 1:
		return sequence_on_direction(myGraph.segment(node_list[0]).sequence, direction_list[0])[start_pos-1:end_pos]

	#Add first node
	found_amr_seq = sequence_on_direction(myGraph.segment(node_list[0]).sequence, direction_list[0])[start_pos-1:]
	#Add middle nodes
	if len(node_list) > 2:
		for i in range(len(node_list)-2):
			seq = sequence_on_direction(myGraph.segment(node_list[i+1]).sequence, direction_list[i+1])
			overlap_size = find_overlap(myGraph.segment(node_list[i]), myGraph.segment(node_list[i+1]))
			seq = seq[overlap_size:]
			found_amr_seq += seq
	#Add last node
	seq = sequence_on_direction(myGraph.segment(node_list[-1]).sequence, direction_list[-1])[:end_pos]
	overlap_size = find_overlap(myGraph.segment(node_list[-2]), myGraph.segment(node_list[-1]))
	seq = seq[overlap_size:]
	found_amr_seq += seq

	return found_amr_seq

def extract_post_sequence_recursively(node, current_seq, current_path, length,
										sequence_list, path_list, node_list):
	"""
	To extract recursively the sequences following the AMR gene sequence and their paths
	Parameters:
		node:		the staring node in the paths
		current_seq:the seq found so far
		current_path: the path found so far
		length: the remained length of the sequences around the AMR gene to be extracted
		sequence_list: the list of sequences following the AMR gene sequence
		path_list:	the list of nodes (and their direction) following the nodes presenting the AMR gene
		node_list:	the list of nodes presenting the AMR gene
	Return:
		modified sequence_list and path_list

	"""
	seq = ""
	found = False
	if length > 0:
		for edge in (node.dovetails_L + node.dovetails_R):
			#the second part of condition ensures that we are finding only nodes that are not presenting the AMR gene itself
			#Logically that shouldn't be the case unless there is a loop in the network:
			# path A -> B -> C represents the AMR; however we have an edge C -> A
			#so when we are looking for the nodes following AMR , A is selected too!!!!!
			if edge.from_segment.name == node.name and edge.to_segment.name != node.name and\
					edge.to_segment.name not in node_list:
				found = True
				new_seq = sequence_on_direction(edge.to_segment.sequence, edge.to_orient)
				#Remove the overlap between nodes' sequences
				overlap_size = find_overlap(node, edge.to_segment)
				new_seq = new_seq[overlap_size:]
				seq = current_seq + new_seq
				path=list.copy(current_path)
				path.append(str(edge.to_segment.name)+str(edge.to_orient))
				if len(new_seq) >= length:
					sequence_list.append(seq[:len(current_seq)+length])
					path_list.append(path)
				else:
					sequence_list, path_list = extract_post_sequence_recursively(edge.to_segment,
					 					seq, path, length - len(new_seq), sequence_list, path_list, node_list)
	if not found:
		sequence_list.append(current_seq)
		path_list.append(current_path)
	return sequence_list, path_list

def extract_pre_sequence_recursively(node, current_seq, current_path,length,
										sequence_list, path_list, node_list):
	"""
	To extract recursively the sequences preceding the AMR gene sequence and their paths
	Parameters:
		node:		the staring node in the paths
		current_seq:the seq found so far
		current_path: the path found so far
		length: the remained length of the sequences around the AMR gene to be extracted
		sequence_list: the list of sequences preceding the AMR gene sequence
		path_list:	the list od nodes (and their direction) preceding the nodes presenting the AMR gene
		node_list:	the list of nodes presenting the AMR gene
	Return:
		modified sequence_list and path_list
	"""
	seq = ""
	found = False
	if length > 0:
		for edge in (node.dovetails_R + node.dovetails_L):
			#the second part ofcondition ensures that we are finding only nodes that are not presenting the AMR gene itself
			#Logically that shouldn't be the case unless there is a loop in the network:
			# path A -> B -> C represents the AMR; however we have an edge C -> A
			#so when we are looking for the nodes preceding AMR , C is selected too!!!!!
			if edge.to_segment.name == node.name and edge.from_segment.name != node.name and\
								edge.from_segment.name not in node_list:
				found = True
				new_seq = sequence_on_direction(edge.from_segment.sequence, edge.from_orient)
				#Remove the overlap between nodes' sequences
				overlap_size = find_overlap(edge.from_segment, node)
				new_seq = new_seq[:len(new_seq)-overlap_size]
				seq =  new_seq + current_seq
				path=list.copy(current_path)
				path.insert(0, str(edge.from_segment.name) + str(edge.from_orient))
				if len(new_seq) >= length:
					sequence_list.append(seq[len(new_seq)-length:])
					path_list.append(path)
				else:
					sequence_list, path_list = extract_pre_sequence_recursively(edge.from_segment,
					 					seq, path, length - len(new_seq), sequence_list, path_list, node_list)
	if not found:
		sequence_list.append(current_seq)
		path_list.append(current_path)
	return sequence_list, path_list

def extract_neighborhood_sequence(gfa_file, length, node_list, direction_list, start_pos, end_pos):
	"""
	To extract all linear sequences +/- length from the start and end of the AMR gene
	Parameters:
		myGraph: The original graph
		length: the length of all sequences around the AMR gene to beextracted
		node_list:	the list of nodes containing the AMR gene (i.e., path)
		direction_list: list of direction of nodes in the path
		start_pos:	the position of the first node from which the AMR gene has started
		end_pos: the position of the last node in which the AMR gene ended
	"""
	myGraph = gfapy.Gfa.from_file(gfa_file)
	prefix_length = postfix_length = length
	pre_sequence_list = []
	post_sequence_list = []
	pre_path_list = []
	post_path_list = []
	#Find the sequence before the AMR gene
	#find the first node in the path and check start_pos
	pre_sequence = ""
	#import pdb; pdb.set_trace()
	if start_pos > 0:
		#attach the beginning of the AMR node (not included in AMR seq) to pre_sequence
		prefix_length -= start_pos
		pre_sequence = sequence_on_direction(myGraph.segment(node_list[0]).sequence, direction_list[0])[:start_pos-1]
		if prefix_length <= 0:
			pre_sequence = pre_sequence[-length:]
			pre_sequence_list.append(pre_sequence)
	#check all edges ended at first node
	if prefix_length > 0:
		pre_sequence_list, pre_path_list = extract_pre_sequence_recursively(
									node = myGraph.segment(node_list[0]),
									current_seq = pre_sequence, current_path = [],
									length = prefix_length, sequence_list = [],
									path_list = [], node_list = node_list)

	#Find the sequence after the AMR gene
	#find the last node in the path and check end_pos
	last_segment = myGraph.segment(node_list[-1])
	post_sequence = ""
	if end_pos > 0:
		#attach the end of the AMR node (not included in AMR seq) to post_sequence
		postfix_length -= (len(last_segment.sequence) - end_pos)
		post_sequence = sequence_on_direction(last_segment.sequence, direction_list[-1])[end_pos:]
		if postfix_length <= 0:
			post_sequence = post_sequence[:length]
			post_sequence_list.append(post_sequence)
	#check all edges started from last_segment
	if postfix_length > 0:
		post_sequence_list, post_path_list = extract_post_sequence_recursively(
									node = last_segment,
									current_seq = post_sequence, current_path = [],
									length = postfix_length, sequence_list = [],
									path_list = [], node_list = node_list)
	#extract found AMR sequence
	found_amr_seq = extract_found_amr(myGraph, node_list, direction_list, start_pos, end_pos)

	#Generating all found sequences
	sequence_list = []
	for pre_seq in pre_sequence_list:
		for post_seq in post_sequence_list:
			sequence_list.append(pre_seq.upper()+ found_amr_seq.lower() +post_seq.upper())

	#add direction to nodes representing AMR gene
	node_dir_list = [node + dir for node, dir in zip(node_list, direction_list)]
	#annotate the path representing AMR gene
	node_dir_list[0]='['+node_dir_list[0]
	node_dir_list[-1]=node_dir_list[-1]+']'
	#Generating all found paths
	path_list = generate_all_paths(pre_path_list, post_path_list, node_dir_list)

	return sequence_list, path_list


def find_amr_related_nodes(amr_file, gfa_file, output_dir, bandage = BANDAGE_PATH):
	"""
	Start with path_nodes=1, if the coverage =100% returen the path!
	if the coverage is less than 100%, increase node_path by one
	and check the coverage
	keep doing that until coverage=100% OR node_path>=MAX_NODE_PATH OR the length
	of the actual path is less than the node_path value
	Then return the coverage and the list of nodes in the path
	"""
	#node_list=set() #first, I thought I collect all nodes from different node_paths but then decided to only keep the last one
	old_nodes = []
	nodes = []
	old_coverage = 0
	coverage = 0
	direction_list = []
	start_pos = 0
	end_pos = 0
	cell_info = ""
	for path_nodes in range(MAX_PATH_NODES):
		output_name=output_dir+'align_result_'+str(path_nodes+1)+'_'+datetime.datetime.now().strftime('%Y-%m-%d_%H-%M-%S')
		if os.path.isfile(output_name+'.tsv'):
			os.remove(output_name+'.tsv')
		command = bandage +' querypaths '+gfa_file+' '+amr_file+' '+output_name + \
			' --pathnodes '+str(path_nodes+1)
		os.system(command)
		print("NOTE: The alignment result from Bandage+BLAST has been stroed in " + output_name+".tsv")
		with open(output_name+".tsv") as tsvfile:
			reader = csv.reader(tsvfile, delimiter='\t')
			#How many rows should I read?????????????????
			#skip the header			
			next(reader)
			for row in reader:
				cell_info = row[1].strip()
				old_coverage = coverage
				coverage = float(re.sub('[%]','',row[3]))
				if int(coverage)==100:
					#row[1].split('(', 1)[0]: to remove (70) from '69-, 2193+ (70)'
					print("NOTE: AMR gene is embeded in the following path: "+cell_info)
					nodes, direction_list, start_pos, end_pos = extract_nodes_in_path(cell_info)
					return nodes, direction_list, start_pos, end_pos, coverage
				else:
					if path_nodes>0:
						old_nodes=list.copy(nodes)
					#nodes, direction_list, start_pos, end_pos = extract_nodes_in_path(row[1].rsplit('(')[0])
					nodes, direction_list, start_pos, end_pos = extract_nodes_in_path(cell_info)
					# if (coverage > old_coverage):
						# for node in nodes:
						# 	node_list.add(node)
					#if nodes and old_nodes are identical; meaning increasing node_paths
					#hasn't changed anything
					if collections.Counter(old_nodes)==collections.Counter(nodes):
						print("AMR gene is embeded in the following path: "+cell_info)
						return nodes, direction_list, start_pos, end_pos, coverage,
						#return node_list, old_coverage
				break #only interested in the first row for now
	print("AMR gene is embeded in the following path: "+cell_info)
	return nodes, direction_list, start_pos, end_pos, coverage

def neighborhood_sequence_extraction(amr_file, gfa_file, length, output_dir, bandage = BANDAGE_PATH):
	"""
	"""
	node_list, direction_list, start_pos, end_pos, coverage = \
		find_amr_related_nodes(amr_file, gfa_file, output_dir, bandage)
	#Extract the sequenc of AMR neighborhood
	sequence_list, path_list = extract_neighborhood_sequence(gfa_file, length,
									node_list, direction_list, start_pos, end_pos)
	seq_file = output_dir+'assembly_neighborhood_sequences_'+str(length)+'_'+\
		datetime.datetime.now().strftime('%Y-%m-%d_%H-%M')+'.txt'
	write_sequences_to_file(sequence_list, path_list, node_list, length, seq_file)
	print("NOTE: The list of neighborhood sequences (extracted from assembly graph)\
	 		has been stroed in " + seq_file)
	return seq_file

def neighborhood_graph_extraction(amr_file, gfa_file, distance, output_dir, bandage = BANDAGE_PATH):
	"""
	"""
	#Retrieve the AMR gene sequence
	amr_seq = retrieve_AMR(amr_file)

	node_list, direction_list, start_pos, end_pos, coverage = \
		find_amr_related_nodes(amr_file, gfa_file, output_dir, bandage)

	#Find neighbourhood of found nodes AND create a new graph contating the target subgraph

	#Find the list of segments and edges
	myGraph, segment_list, edge_list =extract_subgraph(node_list, distance, gfa_file)
	#Generate a graph with the lists of segments and edges
	subGraph = create_graph(myGraph, segment_list, edge_list)
	#write the subgraph in a GFA file
	subgraph_file = output_dir+'neighborhood_subgraph_'+str(distance)+'_'+\
		datetime.datetime.now().strftime('%Y-%m-%d_%H-%M')+'.gfa'
	subGraph.to_file(subgraph_file)
	print("NOTE: The extracted subgraph has been stored in " + subgraph_file)
	return subgraph_file

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
		print('please enter the path for the AMR gene file')
		sys.exit()

	if not args.gfa:
		print('please enter the path for the assembly file')
		sys.exit()

	neighborhood_graph_extraction(arg.amr,args.gfa, args.distance, args.output_dir)
	neighborhood_sequence_extraction(arg.amr,args.gfa, args.length, args.output_dir)
