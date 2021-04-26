
# from extract_neighborhood import extract_amr_neighborhood_in_ref_genome
# from utils import retrieve_AMR
#
# file_path = '/media/Data/PostDoc/Dalhousie/Work/Test2/Experiments/CAMI_M_2/AMR_info/sequences/aadA23.fasta'
# contig_file = '/media/Data/PostDoc/Dalhousie/Work/Test2/Experiments/CAMI_M_2/spade_output/contigs.fasta'
#
# amr_seq, amr_name = retrieve_AMR(file_path)
# import pdb; pdb.set_trace()
# seq_list, contig_name_list = extract_amr_neighborhood_in_ref_genome(amr_seq, contig_file, 1000, 95)
#
# """
# he default blastn command runs megablastn with a default word size of 28,
# whereas '-task blastn' uses a default word size of 11.
# In finding AMRs in ref as well as contigs we use -task blastn while it's not used in the graph
# """

# import gfapy
# import copy
#
# graph_file = '/media/Data/PostDoc/Dalhousie/Work/Test2/test.gfa'
# bcalm_graph = gfapy.Gfa.from_file(graph_file)
# subGraph = gfapy.Gfa(version='gfa1')
# node_list = ['1', '2', '3', '4', '5']
# import pdb; pdb.set_trace()
# for segment_name in node_list:
#     print('segment_name = '+segment_name)
#     seg = bcalm_graph.segment(segment_name)
#     subGraph.append(str(seg))
#     for edge in (seg.dovetails_R + seg.dovetails_L):
#         if edge not in subGraph.edges and\
#         ((edge.to_segment.name == seg.name and edge.from_segment.name in node_list) or\
#         (edge.from_segment.name == seg.name and edge.to_segment.name in node_list)):
#             print('edge: '+edge.from_segment.name+','+edge.to_segment.name)
#             #if edge.is_connected():
#             #    edge.disconnect()
#             subGraph.append(str(edge))
#     #seg.disconnect()

#############################mapping ids############################
# import os
# import sys
# import csv
# import pandas as pd
# import datetime
# import logging
# import collections
# import gzip
# import gfapy
#
# import params
# from extract_neighborhood import neighborhood_sequence_extraction
# from utils import extract_files, restricted_amr_name_from_modified_name,\
#         read_ref_annotations_from_db, amr_name_from_comment, retrieve_AMR
# from full_pipeline import neighborhood_annotation, check_coverage_consistency_remove_rest_seq,\
#         evaluate_sequences_up_down_based_on_coverage
#
# SEQ_NAME_PREFIX = 'ng_sequences_'
# ANNOTATION_DIR = 'annotations'
# EVAL_DIR = 'evaluation'
# WORK_DIR = '111_k31'
# SEARCH_DIR = WORK_DIR+'_r10_search_oh0'
#
# def find_mapping(map_file='mapping.fa'):
#     """
#     read map_file and find the mapping between node ids in spacegraphcats and
#     original bcalm graph
#     map[spacegraphcats_id] = bcalm_id
#     """
#     #map_ids = []
#     map_list = collections.defaultdict(lambda: -1)
#     counter = 0
#     with open(map_file) as myfile:
#         for line in myfile:
#             if line.startswith('>'):
#                 sgc_id = int(line[1:-1])
#                 map_list[sgc_id] = counter
#                 counter+=1
#     return map_list
#
# def extract_neighborhood_graph(node_list, bcalm_graph):
#     """
#     """
#     subGraph = gfapy.Gfa(version='gfa1')
#     for segment_name in node_list:
#         #print('segment_name = '+segment_name)
#         seg = bcalm_graph.segment(segment_name)
#         if seg not in subGraph.segments:
#             subGraph.append(str(seg))
#         for edge in (seg.dovetails_R + seg.dovetails_L):
#             # if edge not in subGraph.edges and\
#             # ((edge.to_segment.name == seg.name and edge.from_segment.name in node_list) or\
#             # (edge.from_segment.name == seg.name and edge.to_segment.name in node_list)):
#             if edge not in subGraph.edges:
#                 subGraph.append(str(edge))
#                 if edge.to_segment.name == seg.name and edge.from_segment not in subGraph.segments:
#                     subGraph.append(str(edge.from_segment))
#                 elif edge.from_segment.name == seg.name and edge.to_segment not in subGraph.segments:
#                     subGraph.append(str(edge.to_segment))
#                 #print('edge: '+edge.from_segment.name+','+edge.to_segment.name)
#     return subGraph
#
# def find_corresponding_amr_file(amr_name, amr_files):
#     """
#     """
#     for amr_file in amr_files:
#         if os.path.splitext(os.path.basename(amr_file))[0]==amr_name:
#             return amr_file
#     print("Error: No amr_file for amr "+amr_name+" was found!")
#     sys.exit()
#
# def main():
#     logging.info("Startting the pipeline for testing spacegraphcats ...")
#     working_dir = params.output_dir+WORK_DIR
#     #find the mapping
#     print('Load mapping ...')
#     map_file = working_dir+'/mapping.fa'
#     map_list = find_mapping(map_file)
#     # extract the list of files
#     searching_dir = params.output_dir + SEARCH_DIR
#     id_files = [os.path.join(searching_dir, f) for f in os.listdir(searching_dir) \
#             if os.path.isfile(os.path.join(searching_dir, f)) and f.endswith('cdbg_ids.txt.gz')]
#     id_dir = params.output_dir+'id_lists_mapped/'
#     if not os.path.exists(id_dir):
#         os.makedirs(id_dir)
#     print('Extarcting the list of nodes for each AMR')
#     graph_node_lists = collections.defaultdict(list)
#     for id_file in id_files:
#         amr_name = os.path.basename(id_file).split('.fasta')[0]
#         with open(id_dir+amr_name+'_id_list.txt', 'w') as ifd:
#             with gzip.open(id_file, 'rt') as fd:
#                 for line in fd:
#                     sgc_id = int(line[:-1])
#                     #convert spc_id to the node_id in the original bcalm graph
#                     bcalm_id = map_list[sgc_id]
#                     graph_node_lists[amr_name].append(str(bcalm_id))
#                     ifd.write(str(bcalm_id)+', ')
#
#     frontier_files = [os.path.join(searching_dir, f) for f in os.listdir(searching_dir) \
#             if os.path.isfile(os.path.join(searching_dir, f)) and f.endswith('frontier.txt.gz')]
#     frontier_dir = params.output_dir+'frontier_lists_mapped/'
#     if not os.path.exists(frontier_dir):
#         os.makedirs(frontier_dir)
#     print('Extarcting the list of nodes for each AMR')
#     graph_node_lists = collections.defaultdict(list)
#     for frontier_file in frontier_files:
#         amr_name = os.path.basename(frontier_file).split('.fasta')[0]
#         with open(frontier_dir+amr_name+'_frontier_list.txt', 'w') as ifd:
#             with gzip.open(frontier_file, 'rt') as fd:
#                 for line in fd:
#                     sgc_id = int(line[:-1])
#                     #convert spc_id to the node_id in the original bcalm graph
#                     bcalm_id = map_list[sgc_id]
#                     ifd.write(str(bcalm_id)+', ')
#
#
#
# if __name__=="__main__":
#     text = 'This code is used to test spacegraphcats'
#     main()


from dna_features_viewer import GraphicFeature, GraphicRecord
import matplotlib.pyplot as plt
import numpy as np
from PIL import Image

seq_info = [{'name':'tnpR', 'start_pos': 557, 'end_pos': 910, 'coverage':47.84},
            {'name':'TEM-1', 'start_pos':1093, 'end_pos':1953, 'coverage':37.44},
            {'name':'APH(6)-Id', 'start_pos':2474, 'end_pos':2064, 'coverage':237.34},
            {'name':'dfrA14', 'start_pos':3017, 'end_pos':2544, 'coverage':12.79}]
features = []
for gene_info in seq_info:
	features.append(GraphicFeature(start = gene_info['start_pos'], end=gene_info['end_pos'],
		strand=+1, color='#ffd700', label = gene_info['name']+'('+str(gene_info['coverage'])+')'))
record = GraphicRecord(sequence_length=3053, features=features)
ax, _  = record.plot(figure_width=10)
ax.figure.savefig('temp.jpg', bbox_inches='tight')
img = Image.open('temp.jpg')
