import sys
import os
import csv
from csv import DictReader

family_file = '/media/Data/PostDoc/Dalhousie/Work/Test2/aro_index.tsv'
amr_list = '/media/Data/PostDoc/Dalhousie/Work/Test2/nucleotide_fasta_protein_homolog_model.fasta'
output_file = '/media/Data/PostDoc/Dalhousie/Work/Test2/nucleotide_fasta_protein_homolog_model_without_efflux.fasta'
efflux_families = ['major facilitator superfamily (MFS) antibiotic efflux pump',
                    'small multidrug resistance (SMR) antibiotic efflux pump',
                    'multidrug and toxic compound extrusion (MATE) transporter',
                    'ATP-binding cassette (ABC) antibiotic efflux pump',
                    'resistance-nodulation-cell division (RND) antibiotic efflux pump']
gene_file = '/media/Data/PostDoc/Dalhousie/Work/Test2/AMR_efflux.txt'
gene_list = []
with open(family_file, 'r') as myfile:
    myreader = DictReader(myfile, delimiter='\t')
    for row in myreader:
        amr_family_list = [x.strip() for x in row['AMR Gene Family'].split(';')]
        for amr_family in amr_family_list:
            if amr_family in efflux_families:
                gene = row['Model Name'].strip().replace(' ','_')
                #gene_name = ''.join(e for e in gene if e.isalpha() or e.isnumeric() or e=='_' or e=='-')
                gene_list.append(gene)
                break

with open(gene_file,'w') as fd:
    for gene in gene_list:
        fd.write(gene+'\n')

myfile = open(output_file, 'w')
with open(amr_list) as fp:
    for line in fp:
        if line.startswith('>'):
            write_next = False
            amr_name = line.split('[')[0].split('|')[-1].strip().replace(' ','_')
            #amr_name_processed = ''.join(e for e in amr_name if e.isalpha() or e.isnumeric() or e=='_' or e=='-')
            if amr_name not in gene_list:
                myfile.write(line)
                write_next = True
        elif write_next:
            myfile.write(line)
myfile.close()
