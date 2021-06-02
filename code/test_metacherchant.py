
import os
import sys
import csv
import pandas as pd
import datetime
import logging

import params
from extract_neighborhood import neighborhood_sequence_extraction
from utils import extract_files, restricted_amr_name_from_modified_name,\
        read_ref_annotations_from_db, amr_name_from_comment, initialize_logger
from full_pipeline import neighborhood_annotation, check_coverage_consistency_remove_rest_seq,\
        evaluate_sequences_up_down_based_on_coverage

SEQ_NAME_PREFIX = 'ng_sequences_'
ANNOTATION_DIR = 'annotations'
EVAL_DIR = 'evaluation'

def extract_amr_info_from_file(amr_info_file):
    """
    """
    amr_names = []
    amr_seqs = []
    with open(amr_info_file) as fd:
        for line in fd:
            if line.startswith('>'):
                amr_comment = line[1:-1]
                continue
            amr_name = amr_name_from_comment(amr_comment)
            amr_names.append(amr_name)
            amr_seqs.append(line[:-1])
    return amr_names, amr_seqs

def find_corresponding_amr_file(amr_name, amr_files):
    """
    """
    for amr_file in amr_files:
        if os.path.splitext(os.path.basename(amr_file))[0]==amr_name:
            return amr_file
    logging.error("Error: No amr_file for amr "+amr_name+" was found!")
    sys.exit()

def main():
    """
    """
    amr_info_file = params.output_dir+'AMR_seqs.fasta'
    amr_names, amr_seqs = extract_amr_info_from_file(amr_info_file)
    ref_amr_files = extract_files(params.amr_files, 'please provide the address of the AMR gene(s)')
    if params.ref_genomes_available:
        df = pd.read_csv(params.ref_ng_annotations_file, skipinitialspace=True,  keep_default_na=False)
        amr_groups = df.groupby('target_amr')
    sequence_dir = params.output_dir+'sequences/'
    if not os.path.exists(sequence_dir):
        os.makedirs(sequence_dir)
    evaluation_dir = params.output_dir+EVAL_DIR+'/'+EVAL_DIR+'_'+str(params.seq_length)+'/'
    if not os.path.exists(evaluation_dir):
        os.makedirs(evaluation_dir)
    summary_file = evaluation_dir+'summaryMetrics_up_down_'+\
        datetime.datetime.now().strftime('%Y-%m-%d_%H-%M')+'.csv'
    with open(summary_file,'a') as fd:
        writer = csv.writer(fd)
        writer.writerow(['AMR', 'Unique_TP#', 'FP#', 'Unique_True#', 'found#','sensitivity', 'precision'])

    average_precision = 0
    average_sensitivity = 0
    no_align_file = open(params.output_dir+'no_align_found.txt', 'w')
    for i, amr_name in enumerate(amr_names):
        logging.info('Processing '+amr_name+' ...')
        restricted_amr_name = restricted_amr_name_from_modified_name(amr_name)
        amr_file = find_corresponding_amr_file(restricted_amr_name, ref_amr_files)
        #sequence extraction
        logging.info('Neighborhood Extraction ...')
        gfa_file = params.output_dir+'output/'+str(i+163)+'/graph.gfa'
        logging.info('gfa_file: '+gfa_file)
        seq_file, path_info_file = neighborhood_sequence_extraction(gfa_file,
                params.seq_length, sequence_dir, params.BANDAGE_PATH,
                params.amr_identity_threshold, SEQ_NAME_PREFIX+str(i+163),
                params.path_node_threshold , params.path_seq_len_percent_threshod,
                params.max_kmer_size, params.assembler, (amr_file,''))
        if seq_file=='':
            no_align_file.write(amr_name+'\n')
            sensitivity = 0
            precision = 0
            with open(summary_file,'a') as fd:
                writer = csv.writer(fd)
                writer.writerow([amr_name, 0, 0, 0, -1,0, 0])
        else:
            #Annotation
            logging.info('Annoation ...')
            if params.ref_genomes_available:
                ref_up_info_list, ref_amr_info_list, ref_down_info_list =\
                    read_ref_annotations_from_db(amr_groups, amr_name)
            all_seq_info_list, annotation_file =neighborhood_annotation(amr_name, seq_file,
                    path_info_file, params.seq_length, [], [], params.output_dir,
                    params.PROKKA_COMMAND_PREFIX,params.use_RGI,
                    params.RGI_include_loose, '_'+str(i+163)+restricted_amr_name,
                    params.amr_identity_threshold, False)
            #Evaluation
            logging.info('Evaluation ...')
            sensitivity, precision = evaluate_sequences_up_down_based_on_coverage(
                    amr_name, annotation_file, summary_file,ref_up_info_list,
                    ref_down_info_list, ref_amr_info_list,params.assembler)
        # annotate_dir = params.output_dir+ANNOTATION_DIR+'/'+ANNOTATION_DIR+'_'+\
		# 	str(params.seq_length)+'/annotation_'+restricted_amr_name+'_'+str(params.seq_length)
		# coverage_annotation, remained_seqs = check_coverage_consistency_remove_rest_seq(\
		# 				all_seq_info_list, params.output_dir,params.coverage_thr,
        #                 restricted_amr_name, annotate_dir+'/')
        # #Evaluation
        # print('Evaluation ...')
        # sensitivity, precision = evaluate_sequences_up_down_based_on_coverage(
        #         amr_name, coverage_annotation, summary_file,
        #         ref_up_info_list, ref_down_info_list, ref_amr_info_list)
        average_precision+=precision
        average_sensitivity+=sensitivity
        logging.info('For "'+amr_name+'": sensitivity= '+str(sensitivity)+' precision = '+ str(precision))

    average_precision = average_precision/len(amr_names)
    average_sensitivity = average_sensitivity/len(amr_names)
    logging.info("average_precision = "+str(average_precision))
    logging.info("average_sensitivity = "+str(average_sensitivity))
    no_align_file.close()
    # else:
    #     print('please enter the path for amr sequence files')
    #     import pdb; pdb.set_trace()
    #     sys.exit()


if __name__=="__main__":
    text = 'This code is used to test metacherchant'
    log_name = 'logger_'+datetime.datetime.now().strftime('%Y-%m-%d_%H-%M')+'.log'
    initialize_logger(params.main_dir, log_name)
    # parser = argparse.ArgumentParser(description=text)
    # parser.add_argument('--amr', type=str, default='',
    #     help='the path of the directory containing the amr files')
    # parser.add_argument('--seq', type=str, default = '',
    #     help = 'the path of the fasta file containing all AMR sequences')
    # args = parser.parse_args()
    main()
