"""
File:		generate_upset_csv.py
Aythor:		Somayeh Kafaie
Date:		March 2021


Purpose:	To extract the list of AMRs that their neighborhood has been
            detected successfully across different settings (and generate upset plot)

To run:
	generate_upset_csv.py
"""

import csv
from csv import DictReader
import collections
from upsetplot import plot
from upsetplot import from_memberships
from matplotlib import pyplot
import pandas as pd

from utils import extract_files, retrieve_AMR

amr_files = 'AMR_info/sequences/'
metaspade_file ='comparison/cami_h_1_metaspades_seqLen1000_Cov30.csv'
#bcalm_file =
contig_file ='comparison/cami_h_1_contigs_seqLen1000.csv'
file_list = [metaspade_file, contig_file]
sensitivity_thr_list  = [0.5, 1]

def read_sensitivity_info_from_file(summary_file, amr_ng_presence, sensitivity_thr):
    """
    """
    with open(summary_file, 'r') as myfile:
        myreader = DictReader(myfile, delimiter=',')
        for row in myreader:
            if row['AMR']!='':
                if float(row['sensitivity'].strip()) >=sensitivity_thr:
                    amr_ng_presence[row['AMR']].append(1)
                else:
                    amr_ng_presence[row['AMR']].append(0)
    return amr_ng_presence

def main():
    """
    """
    for sensitivity_thr in sensitivity_thr_list:
        csv_file = 'comparison/AMR_presence_sensitivity_'+str(sensitivity_thr)+'.csv'
        #ref_amr_files = extract_files(amr_files, 'please provide the address of the AMR gene(s)')
        #ref_amr_names = [retrieve_AMR(e)[1] for e in ref_amr_files]
        amr_ng_presence = collections.defaultdict(list)
        for summary_file in file_list:
            amr_ng_presence = read_sensitivity_info_from_file(summary_file,
                                        amr_ng_presence, sensitivity_thr)
        #set_list = []
        columns = ['Reference','MetaSpades','Contigs']
        #write into final csv file
        with open(csv_file, 'w') as fd:
            writer = csv.writer(fd)
            writer.writerow(['AMR_name']+columns)
            for k, v in amr_ng_presence.items():
                #values = ','.join(e for e in v)
                writer.writerow( [k] + [1] + v)
                #set_list.append([columns[i] for i, e in enumerate([1]+v) if e==1 ])
        # test = from_memberships(set_list)
        # import pdb; pdb.set_trace()
        # plot(test)
        # pyplot.show()


if __name__=="__main__":
    text = 'This code is used to extract the list of AMRs that their neighborhood'+\
            ' has been detected successfully across different settings (and generate'+\
            ' upset plot)'
    main()
