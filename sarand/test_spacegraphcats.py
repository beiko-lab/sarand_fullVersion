import os
import sys
import csv
import pandas as pd
import datetime
import logging
import collections
import gzip
import gfapy

import params
from extract_neighborhood import neighborhood_sequence_extraction
from utils import extract_files, restricted_amr_name_from_modified_name,\
        read_ref_annotations_from_db, amr_name_from_comment, retrieve_AMR,\
        initialize_logger
from full_pipeline import neighborhood_annotation, check_coverage_consistency_remove_rest_seq,\
        evaluate_sequences_up_down_based_on_coverage

SEQ_NAME_PREFIX = 'ng_sequences_'
ANNOTATION_DIR = 'annotations'
EVAL_DIR = 'evaluation'
WORK_DIR = '222_k41'
SEARCH_DIR = WORK_DIR+'_r10_search_oh0'

def find_mapping(map_file='mapping.fa'):
    """
    read map_file and find the mapping between node ids in spacegraphcats and
    original bcalm graph
    map[spacegraphcats_id] = bcalm_id
    """
    #map_ids = []
    map_list = collections.defaultdict(lambda: -1)
    counter = 0
    with open(map_file) as myfile:
        for line in myfile:
            if line.startswith('>'):
                sgc_id = int(line[1:-1])
                map_list[sgc_id] = counter
                counter+=1
    return map_list

def extract_neighborhood_graph(node_list, bcalm_graph):
    """
    """
    subGraph = gfapy.Gfa(version='gfa1')
    for segment_name in node_list:
        #print('segment_name = '+segment_name)
        if segment_name == '1067521':
            import pdb; pdb.set_trace()
        seg = bcalm_graph.segment(segment_name)
        if seg not in subGraph.segments:
            subGraph.append(str(seg))
        for edge in (seg.dovetails_R + seg.dovetails_L):
            # if edge not in subGraph.edges and\
            # ((edge.to_segment.name == seg.name and edge.from_segment.name in node_list) or\
            # (edge.from_segment.name == seg.name and edge.to_segment.name in node_list)):
            if edge not in subGraph.edges:
                subGraph.append(str(edge))
                if edge.to_segment.name == seg.name and edge.from_segment not in subGraph.segments:
                    subGraph.append(str(edge.from_segment))
                elif edge.from_segment.name == seg.name and edge.to_segment not in subGraph.segments:
                    subGraph.append(str(edge.to_segment))
                #print('edge: '+edge.from_segment.name+','+edge.to_segment.name)
    #add edges that are between segments that are not in the node_list but in their direct connection
    for segment in subGraph.segments:
        if segment.name not in node_list:
            seg = bcalm_graph.segment(segment.name)
            for edge in seg.dovetails:
                if edge.from_segment in subGraph.segments and\
                    edge.to_segment in subGraph.segments and edge not in subGraph.edges:
                    subGraph.append(str(edge))
    return subGraph

def find_corresponding_amr_file(amr_name, amr_files):
    """
    """
    for amr_file in amr_files:
        if os.path.splitext(os.path.basename(amr_file))[0]==amr_name:
            return amr_file
    logging.error("Error: No amr_file for amr "+amr_name+" was found!")
    sys.exit()

def main():
    logging.info("Startting the pipeline for testing spacegraphcats ...")
    working_dir = params.output_dir+WORK_DIR
    #converting graph to gfa format
    logging.info('converting graph to gfa format ... ')
    command = 'python /media/Data/tools/bcalm_convertToGFA.py '+\
    working_dir+'/bcalm.unitigs.fa '+working_dir+'/'+params.gfa_file+' 31'
    os.system(command)
    logging.info('Loading graph from gfa file ...')
    bcalm_graph = gfapy.Gfa.from_file(working_dir+'/'+params.gfa_file)
    #find the mapping
    logging.info('Load mapping ...')
    map_file = working_dir+'/mapping.fa'
    map_list = find_mapping(map_file)
    # extract the list of files
    searching_dir = params.output_dir + SEARCH_DIR
    id_files = [os.path.join(searching_dir, f) for f in os.listdir(searching_dir) \
            if os.path.isfile(os.path.join(searching_dir, f)) and f.endswith('cdbg_ids.txt.gz')]
    # if compressed_files:
    #     for c_file in compressed_files:
    #         command('gunzip '+c_file)
    #         os.system(command)
    # id_files = [os.path.join(gfiles, f) for f in os.listdir(gfiles) \
    # 		if os.path.isfile(os.path.join(gfiles, f)) and f.endswith('cdbg_ids.txt')]
    #Extract the list of nodes representing the sungraph for each AMR
    logging.info('Extarcting the list of nodes for each AMR')
    graph_node_lists = collections.defaultdict(list)
    for id_file in id_files:
        amr_name = os.path.basename(id_file).split('.fasta')[0]
        with gzip.open(id_file, 'rt') as fd:
            for line in fd:
                sgc_id = int(line[:-1])
                #convert spc_id to the node_id in the original bcalm graph
                bcalm_id = map_list[sgc_id]
                graph_node_lists[amr_name].append(str(bcalm_id))
    #extract the list of amr files
    amr_files = extract_files(params.amr_files, 'please provide the address of AMR files')
    if params.ref_genomes_available:
        df = pd.read_csv(params.ref_ng_annotations_file, skipinitialspace=True,  keep_default_na=False)
        amr_groups = df.groupby('target_amr')
    sequence_dir = params.output_dir+'sequences/'
    if not os.path.exists(sequence_dir):
        os.makedirs(sequence_dir)
    evaluation_dir = params.output_dir+EVAL_DIR+'/'+EVAL_DIR+'_'+str(params.seq_length)+'/'
    if not os.path.exists(evaluation_dir):
        os.makedirs(evaluation_dir)
    summary_file = evaluation_dir+'summaryMetrics_up_down_'+datetime.datetime.now().strftime('%Y-%m-%d_%H-%M')+'.csv'
    with open(summary_file,'a') as fd:
        writer = csv.writer(fd)
        writer.writerow(['AMR', 'Unique_TP#', 'FP#', 'Unique_True#', 'found#','sensitivity', 'precision'])
    graph_dir = params.output_dir+'ng_graphs/'
    if not os.path.exists(graph_dir):
        os.makedirs(graph_dir)
    average_precision = 0
    average_sensitivity = 0
    no_align_file = open(params.output_dir+'no_align_found.txt', 'w')
    #Going over the list of AMRs
    for r_amr_name, node_list in graph_node_lists.items():
        logging.info('Processing '+r_amr_name+' ...')
        #extract graph from the nodes
        ng_graph = extract_neighborhood_graph(node_list, bcalm_graph)
        gfa_file = graph_dir+r_amr_name+'.gfa'
        ng_graph.to_file(gfa_file)
        amr_file = find_corresponding_amr_file(r_amr_name, amr_files)
        amr_seq, amr_name = retrieve_AMR(amr_file)
        #sequence extraction
        logging.info('Neighborhood Extraction ...')
        seq_file, path_info_file = neighborhood_sequence_extraction(gfa_file,
                params.seq_length, sequence_dir, params.BANDAGE_PATH,
                params.amr_identity_threshold, SEQ_NAME_PREFIX+'_'+r_amr_name,
                params.path_node_threshold , params.path_seq_len_percent_threshod,
                params.max_kmer_size, params.assembler, (amr_file,''))
        if seq_file=='':
            no_align_file.write(amr_name+'\n')
            sensitivity = 0
            precision = 0
            with open(summary_file,'a') as fd:
                writer = csv.writer(fd)
                writer.writerow([amr_name, 0, 0, 0, -1, 0, 0])
        else:
            #Annotation
            logging.info('Annoation ...')
            if params.ref_genomes_available:
                ref_up_info_list, ref_amr_info_list, ref_down_info_list =\
                    read_ref_annotations_from_db(amr_groups, amr_name)
            all_seq_info_list, annotation_file =neighborhood_annotation(amr_name, seq_file,
                    path_info_file, params.seq_length, [], [], params.output_dir,
                    params.PROKKA_COMMAND_PREFIX,params.use_RGI,
                    params.RGI_include_loose, '_'+r_amr_name,
                    params.amr_identity_threshold, False)
            #Evaluation
            logging.info('Evaluation ...')
            sensitivity, precision = evaluate_sequences_up_down_based_on_coverage(
                    amr_name, annotation_file, summary_file,ref_up_info_list,
                    ref_down_info_list, ref_amr_info_list,params.assembler)
        average_precision+=precision
        average_sensitivity+=sensitivity
        logging.info('For "'+amr_name+'": sensitivity= '+str(sensitivity)+' precision = '+ str(precision))

    average_precision = average_precision/len(amr_files)
    average_sensitivity = average_sensitivity/len(amr_files)
    logging.info("average_precision = "+str(average_precision))
    logging.info("average_sensitivity = "+str(average_sensitivity))
    no_align_file.close()

if __name__=="__main__":
    text = 'This code is used to test spacegraphcats'
    log_name = 'logger_'+datetime.datetime.now().strftime('%Y-%m-%d_%H-%M')+'.log'
    initialize_logger(params.main_dir, log_name)
    main()
