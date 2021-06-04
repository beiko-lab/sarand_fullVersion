"""
File:		params.py
Aythor:		Somayeh Kafaie
Date:		September 2020

Purpose:	To set the parameters used in the project

NOTE: all parameters can be set in params.py
NOTE: if use_RGI = TRUE, make sure either RGI has been installed system-wide or
	you already are in the environment RGI installed in!
"""
import enum

"""
# The list of valid tasks
# to use them set task in params.py:
# for a single task: insert its value as the single item in task list
# for a range of tasks: insert the value of start and end as the first and the second (last) item in the list
# for all tasks: task = [0]
"""
class Pipeline_tasks(enum.Enum):
	all = 0
	metagenome_creation = 1
	read_simulation = 2
	assembly = 3
	graph_neighborhood = 4
	sequence_neighborhood = 5
	neighborhood_annotation = 6
	neighborhood_evaluation = 7

# Whether to insert the AMR sequence in arandom or pre-defined location of the genome
# in case of amr artificial insertion
class Insertion_type(enum.Enum):
	random = 1
	assigned = 2

class Assembler_name(enum.Enum):
    metaspades = 1
    megahit = 2
    bcalm = 3
    metacherchant = 4
    spacegraphcats = 5

BANDAGE_PATH = 'Bandage'
#BANDAGE_PATH = '/media/Data/tools/Bandage_Ubuntu_dynamic_v0_8_1/Bandage'
ART_PATH = 'art_illumina'
SPADES_PATH = 'spades.py'
PROKKA_COMMAND_PREFIX = 'docker run -v `pwd`:/data staphb/prokka:latest '

CARD_AMR_SEQUENCES ='data/CARD_AMR_seq.fasta'
main_dir = 'test/'
output_dir = main_dir+'output_dir/'

multi_processor = True
core_num = 4
coverage_thr = 30
task = [0]

amr_files = output_dir +'AMR_info/sequences/'
find_amr_genes = False
amr_identity_threshold =  95

ref_genomes_available = True
ref_ng_annotations_file = output_dir + 'AMR_info/ref_neighborhood_annotations.csv'
ref_genome_files = main_dir + 'metagenome_data/'


#simulating reads
read_length =  150
metagenome_file = output_dir +'metagenome.fasta'
#reads =[output_dir + 'sub50_trimmed_ERR1713331_1.fastq',\
# 		output_dir + 'sub50_trimmed_ERR1713331_2.fastq']
#reads = output_dir +'H_S001__insert_180_reads_anonymous.fq.gz'
reads =[output_dir + 'NW016_1_trimmed._R1.fastq',\
 		output_dir + 'NW016_1_trimmed._R2.fastq']

spades_thread_num = 16

assembler = Assembler_name.metaspades
#Setting for assembler
if assembler == Assembler_name.metaspades:
	assembler_output_dir = 'spade_output'
	gfa_file = output_dir + assembler_output_dir + '/assembly_graph_with_scaffolds.gfa'
	contig_file = output_dir + assembler_output_dir + '/contigs.fasta'
	max_kmer_size = 55
elif assembler == Assembler_name.megahit:
	assembler_output_dir = 'megahit_output'
	gfa_file = output_dir + assembler_output_dir + '/k59.gfa'
	#gfa_file = output_dir + assembler_output_dir + '/k119.gfa'
	#contig_file = output_dir + assembler_output_dir + '/final.contigs.fa'
	contig_file = output_dir + assembler_output_dir + '/intermediate_contigs/k59.contigs.fa'
	#max_kmer_size = 59
	max_kmer_size = 119
elif assembler == Assembler_name.bcalm:
	assembler_output_dir = 'bcalm_output'
	gfa_file = output_dir + assembler_output_dir + '/bcalm_graph_55.gfa'
	#doesn't produce any contig file; so just sue the one from meta-pades for evaluation
	contig_file = output_dir + 'spades_output/contigs.fasta'
	max_kmer_size = 54
elif assembler == Assembler_name.metacherchant:
	max_kmer_size = 30
elif assembler == Assembler_name.spacegraphcats:
	gfa_file = 'bcalm_graph.gfa'
	max_kmer_size = 30

spades_error_correction = True

graph_distance = 3
seq_length = 1000

# to be used in sequence extraction and after having path_node_threshold number of
#nodes in our path or already extracting path_seq_len_percent_threshold percent of
#neighborhood sequence, we stop or only traverse the path from the longest neighbor of the current node
path_node_threshold =  10
path_seq_len_percent_threshold = 90

use_RGI =  True
RGI_include_loose = False

ng_seq_files = output_dir + 'sequences_info/sequences_info_'+str(seq_length)+'/sequences/'
ng_path_info_files = output_dir + 'sequences_info/sequences_info_'+str(seq_length)+'/paths_info/'


#insertion of AMR gene in ref_genome
artificial_amr_insertion = False
inserted_amr_file = ''
number_of_insertions = [2, 1]
number_of_copies = [[3, 1],[1]]
insertion_type = Insertion_type.random
insertion_locations = []

genome_amr_files = [output_dir + "SE_FDAARGOS_768_TRU-1_37589.fasta",
                    output_dir + "enter_zy2_TRU1_19098.fasta"]
