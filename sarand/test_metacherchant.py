
import os
import sys
import csv
import pandas as pd
import datetime
import logging

import params
from params import Pipeline_tasks
from extract_neighborhood import neighborhood_sequence_extraction, read_path_info_from_align_file
from utils import extract_files, restricted_amr_name_from_modified_name,\
        read_ref_annotations_from_db, amr_name_from_comment, initialize_logger,\
        validate_task_values
from full_pipeline import neighborhood_annotation, check_coverage_consistency_remove_rest_seq,\
        evaluate_sequences_up_down_based_on_coverage

SEQ_NAME_PREFIX = 'ng_sequences_'
ANNOTATION_DIR = 'annotations'
EVAL_DIR = 'evaluation'

def find_annotation_file(annotation_file, out_dir,amr_name, seq_len):
    """
    """
    if annotation_file!='':
        return annotation_file
    annotate_dir = os.path.join(out_dir, ANNOTATION_DIR, ANNOTATION_DIR+'_'+str(seq_len),
                                    'annotation_'+amr_name+'_'+str(seq_len))
    for f in os.listdir(annotate_dir):
        if f=='trimmed_annotation_info_'+amr_name+'.csv':
            return os.path.join(annotate_dir, f)
    return ''

def find_sequence_file(seq_file, sequence_dir, amr_name, seq_len):
    """
    """
    #seq_file has already been found
    if seq_file!='':
        return seq_file
    seq_dir = os.path.join(sequence_dir, 'sequences')
    for f in os.listdir(seq_dir):
        if f.startswith(SEQ_NAME_PREFIX+amr_name+'_'+str(seq_len)) and\
            os.path.isfile(os.path.join(seq_dir, f)) and\
            os.path.splitext(f)[1]=='.txt':
            return os.path.join(seq_dir, f)
    return ''

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
    task_list = validate_task_values(params.task)
    amr_info_file = os.path.join(params.main_dir, 'AMR_seqs_full.fasta')
    amr_names, amr_seqs = extract_amr_info_from_file(amr_info_file)
    ref_amr_files = extract_files(params.amr_files, 'please provide the address of the AMR gene(s)')
    if params.ref_genomes_available:
        df = pd.read_csv(params.ref_ng_annotations_file, skipinitialspace=True,  keep_default_na=False)
        amr_groups = df.groupby('target_amr')
    sequence_dir = os.path.join(params.output_dir, 'sequences')
    if not os.path.exists(sequence_dir):
        os.makedirs(sequence_dir)
    ##to make use of previously generated alignment files
    align_files = []
    align_dir = os.path.join(sequence_dir, 'alignment_files')
    if os.path.exists(align_dir):
        align_files = [os.path.join(align_dir, f) for f in os.listdir(align_dir) \
				if os.path.isfile(os.path.join(align_dir, f)) and os.path.splitext(f)[1]=='.tsv']
    evaluation_dir = os.path.join(params.output_dir, EVAL_DIR, EVAL_DIR+'_'+str(params.seq_length))
    if not os.path.exists(evaluation_dir):
        os.makedirs(evaluation_dir)
    summary_file = os.path.join(evaluation_dir, 'summaryMetrics_up_down_metacherchant_'+\
        datetime.datetime.now().strftime('%Y-%m-%d_%H-%M')+'.csv')
    with open(summary_file,'a') as fd:
        writer = csv.writer(fd)
        writer.writerow(['AMR', 'Unique_TP#', 'FP#', 'Unique_True#', 'found#','sensitivity', 'precision', 'not_found'])

    average_precision = 0
    average_sensitivity = 0
    no_align_file = open(os.path.join(params.output_dir, 'no_align_found.txt'), 'w')
    for amr_name in amr_names:
        logging.info('Processing '+amr_name+' ...')
        restricted_amr_name = restricted_amr_name_from_modified_name(amr_name)
        #amr_name1 = amr_name.replace(";",'SS')
        #restricted_amr_name1 = ''.join(e for e in amr_name1 if e.isalpha() or e.isnumeric() or e=='_' or e=='-')
        seq_file = ''
        annotation_file = ''
        # sequence extraction
        if Pipeline_tasks.sequence_neighborhood.value in task_list:
            logging.info('Neighborhood Extraction ...')
            amr_file = find_corresponding_amr_file(restricted_amr_name, ref_amr_files)
            found = False
            for myfile in align_files:
                if os.path.basename(myfile).startswith('ng_sequences_'+restricted_amr_name+'_align') or\
                os.path.basename(myfile).startswith('ng_sequences_'+restricted_amr_name+restricted_amr_name+'_align'):
                    align_file = myfile
                    found, paths_info = read_path_info_from_align_file(align_file, 95)
                    break
            if found:
                amr_pair = (amr_file,paths_info)
            else:
                amr_pair = (amr_file, '')
            #gfa_file = params.output_dir+'output/'+str(i+1)+'/graph.gfa'
            gfa_file = os.path.join(params.main_dir, 'output', restricted_amr_name, 'graph.gfa')
            command = "sed -e 's/\tCL:z:GREEN//g' -i "+gfa_file
            print(command)
            os.system(command)
            logging.info('gfa_file: '+gfa_file)
            seq_file, path_info_file = neighborhood_sequence_extraction(gfa_file,
                    params.seq_length, sequence_dir, params.BANDAGE_PATH,
                    params.amr_identity_threshold, SEQ_NAME_PREFIX,
                    params.path_node_threshold , params.path_seq_len_percent_threshold,
                    params.max_kmer_size, params.ng_extraction_time_out,
                    params.assembler, amr_pair)

        #Annotation
        if Pipeline_tasks.neighborhood_annotation.value in task_list:
            logging.info('Annoation ...')
            seq_file = find_sequence_file(seq_file, sequence_dir, restricted_amr_name, params.seq_length)
            if seq_file=='':
                print('no sequence file is available for '+amr_name)
                no_align_file.write(amr_name+'\n')
                sensitivity = 0
                precision = 0
                with open(summary_file,'a') as fd:
                    writer = csv.writer(fd)
                    #writer.writerow([amr_name, 0, 0, 0, -1,0, 0])
                    writer.writerow([amr_name, 0, 0, 0, -1,0, 0, True])
                continue
            _, annotation_file =neighborhood_annotation(amr_name, seq_file,
                    path_info_file, params.seq_length, [], [], params.output_dir,
                    params.PROKKA_COMMAND_PREFIX,params.use_RGI,
                    params.RGI_include_loose, '_'+restricted_amr_name,
                    False)

        #Evaluation
        if Pipeline_tasks.neighborhood_evaluation.value in task_list and\
            params.ref_genomes_available:
            logging.info('Evaluation ...')
            annotation_file = find_annotation_file(annotation_file, params.output_dir,
                restricted_amr_name, params.seq_length)
            if annotation_file =='':
                print('no annotation file is available for '+amr_name)
            ref_up_info_list, ref_amr_info_list, ref_down_info_list =\
                read_ref_annotations_from_db(amr_groups, amr_name)
            sensitivity, precision = evaluate_sequences_up_down_based_on_coverage(
                    amr_name, annotation_file, summary_file,ref_up_info_list,
                    ref_down_info_list, ref_amr_info_list,params.assembler)
            average_precision+=precision
            average_sensitivity+=sensitivity
            logging.info('For "'+amr_name+'": sensitivity= '+str(sensitivity)+' precision = '+ str(precision))

    average_precision = average_precision/len(amr_names)
    average_sensitivity = average_sensitivity/len(amr_names)
    logging.info("average_precision = "+str(average_precision))
    logging.info("average_sensitivity = "+str(average_sensitivity))
    no_align_file.close()


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
