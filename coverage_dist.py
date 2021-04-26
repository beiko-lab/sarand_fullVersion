import os
import csv
from csv import DictReader
import pandas as pd
import collections
import seaborn as sns
import matplotlib.pyplot as plt

from utils import retreive_original_amr_name, extract_up_down_from_csv_file,\
            similar_seq_annotation_already_exist, extract_files, \
            read_info_from_overlap_ref_files, read_ref_annotations_from_db,\
            extract_name_from_file_name, retrieve_AMR
from full_pipeline import exists_in_gene_list


ANNOTATION_DIR = 'annotations'
AMR_DIR_NAME = 'AMR_info/'
AMR_OVERLAP_FILE = 'overlaps.txt'


def evaluate_sequences_gene_based_on_coverage(amr_name, annotation_file,
                        ref_up_info_list, ref_down_info_list, ref_amr_info_list,
                        fp_file_name, tp_file_name):
    """
    """

    fp_file = open(fp_file_name, 'a')
    fp_writer = csv.writer(fp_file)
    tp_file = open(tp_file_name, 'a')
    tp_writer = csv.writer(tp_file)
    original_amr_name = retreive_original_amr_name(amr_name)
    if annotation_file=="":
        return
    #Read  annotation_file csv file and store info
    seq_info = []
    #seq_info_list = []
    up_info_list = []
    down_info_list = []
    amr_info_list = []
    sequence = ''
    amr_coverage = 0
    with open(annotation_file, 'r') as myfile:
        myreader = DictReader(myfile)
        old_seq = ''
        for row in myreader:
            if row['seq_name'].startswith('extracted'):
                if row['coverage']=='':
                    print("ERROR: Coverage information are not available for "+amr_name)
                    import pdb; pdb.set_trace()
                    sys.exit()
                #else:
                gene_info = {'seq_name':row['seq_name'], 'seq_value':row['seq_value'],
                            'gene':row['gene'], 'length':row['length'],
                            'start_pos':int(row['start_pos']),'end_pos':int(row['end_pos']),
                            'target_amr':row['target_amr'], 'coverage':float(row['coverage'])}
                cur_seq = row['seq_name']
                if cur_seq!=old_seq:
                    if (seq_info):
                        amr_found, up_info, down_info, amr_info = extract_up_down_from_csv_file(seq_info)
                        if amr_found:
                            amr_info_list.append(amr_info)
                            amr_coverage = float(amr_info['coverage'])
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
            amr_coverage = float(amr_info['coverage'])
            if up_info and not similar_seq_annotation_already_exist(up_info, up_info_list):
                up_info_list.append(up_info)
            if down_info and not similar_seq_annotation_already_exist(down_info, down_info_list):
                down_info_list.append(down_info)
    #import pdb; pdb.set_trace()
    for seq_info in up_info_list:
        if seq_info:
            for gene_info in seq_info:
                if exists_in_gene_list(ref_up_info_list, gene_info):
                    diff = abs(amr_coverage-gene_info['coverage'])
                    tp_writer.writerow([original_amr_name, gene_info['gene'], amr_coverage,
                                gene_info['coverage'], diff, diff/amr_coverage,
                                diff/(amr_coverage+gene_info['coverage'])])
                else:
                    diff = abs(amr_coverage-gene_info['coverage'])
                    fp_writer.writerow([original_amr_name, gene_info['gene'], amr_coverage,
                                gene_info['coverage'], diff, diff/amr_coverage,
                                diff/(amr_coverage+gene_info['coverage'])])
    for seq_info in down_info_list:
        if seq_info:
            for gene_info in seq_info:
                if exists_in_gene_list(ref_down_info_list, gene_info):
                    diff = abs(amr_coverage-gene_info['coverage'])
                    tp_writer.writerow([original_amr_name, gene_info['gene'], amr_coverage,
                                gene_info['coverage'], diff, diff/amr_coverage,
                                diff/(amr_coverage+gene_info['coverage'])])
                else:
                    diff = abs(amr_coverage-gene_info['coverage'])
                    fp_writer.writerow([original_amr_name, gene_info['gene'], amr_coverage,
                                gene_info['coverage'], diff, diff/amr_coverage,
                                diff/(amr_coverage+gene_info['coverage'])])
    fp_file.close()
    tp_file.close()

def main(params):
    """
    """
    fp_file_name = params.output_dir+'results_metaspades/false_positive_dist.csv'
    tp_file_name = params.output_dir+'results_metaspades/true_positive_dist.csv'
    with open(fp_file_name, 'a') as fp_file:
        writer_fp = csv.writer(fp_file)
        writer_fp.writerow(['target_amr', 'gene', 'amr_coverage', 'gene_coverage',
                'diff_coverage', 'diff_DIV_amr', 'diff_DIV_sum'])
    with open(tp_file_name, 'a') as tp_file:
        writer_tp = csv.writer(tp_file)
        writer_tp.writerow(['target_amr', 'gene', 'amr_coverage', 'gene_coverage',
                'diff_coverage', 'diff_DIV_amr', 'diff_DIV_sum'])
    #to find the annotation of ref genomes for all AMRs
    df = pd.read_csv(params.ref_ng_annotations_file, skipinitialspace=True,  keep_default_na=False)
    amr_groups = df.groupby('target_amr')
    #find the list of amr files
    overlap_file_name = params.output_dir+'results_metaspades/'+AMR_DIR_NAME+AMR_OVERLAP_FILE
    if os.path.isfile(overlap_file_name):
        ref_amr_files = extract_files(params.amr_files, 'please provide the address of the AMR gene(s)')
        unique_amr_files, not_found_amr_names, amr_count = read_info_from_overlap_ref_files(
					overlap_file_name, ref_amr_files)

    for i, amr_file in enumerate(unique_amr_files):
        restricted_amr_name = extract_name_from_file_name(amr_file)
        _, amr_name = retrieve_AMR(amr_file)
        #remove some extracted sequences based on coverage consistency
        annotate_dir = params.output_dir+'results_metaspades/'+ANNOTATION_DIR+'/'+ANNOTATION_DIR+'_'+\
			str(params.seq_length)+'/annotation_'+restricted_amr_name+'_'+str(params.seq_length)
        annotation_file = annotate_dir+'/trimmed_annotation_info_'+restricted_amr_name+'.csv'
        ref_up_info_list, ref_amr_info_list, ref_down_info_list =\
                read_ref_annotations_from_db(amr_groups, amr_name)

        evaluate_sequences_gene_based_on_coverage(amr_name, annotation_file,
                        ref_up_info_list, ref_down_info_list, ref_amr_info_list,
                        fp_file_name, tp_file_name)

    for col in ['diff_coverage','diff_DIV_amr','diff_DIV_sum']:
        output_file = params.output_dir+'results_metaspades/'+'TP_2_2_2_'+col+'.png'
        data = pd.read_csv(tp_file_name, skipinitialspace=True, usecols=[col])
        sns.displot(data, kde=True, rug=False)
        plt.title('True positives (2_2_2)')
        plt.savefig(output_file)
        # g.set_titles('True positivess (2_2_2)')
        # g.savefig(output_file)
    for col in ['diff_coverage','diff_DIV_amr','diff_DIV_sum']:
        output_file = params.output_dir+'results_metaspades/'+'FP_2_2_2_'+col+'.png'
        data = pd.read_csv(fp_file_name, skipinitialspace=True, usecols=[col])
        sns.displot(data, kde=True, rug=False)
        plt.title('False positives (2_2_2)')
        plt.savefig(output_file)
        #g.set_titles('False positives (2_2_2)')
        # fig = g.get_figure()
        # fig.savefig(output_file)
        #g.savefig(output_file)
#data = pd.read_csv(fp_file_name, skipinitialspace=True, usecols=['diff_coverage',
#'diff_DIV_amr','diff_DIV_sum'])

if __name__=="__main__":
	import params
	main(params)
