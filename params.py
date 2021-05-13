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

#from full_pipeline import Pipeline_tasks
#from full_pipeline import Insertion_type

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

BANDAGE_PATH = '/media/Data/tools/Bandage_Ubuntu_dynamic_v0_8_1/Bandage'
ART_PATH = '/media/Data/tools/art_src_MountRainier_Linux/art_illumina'
SPADES_PATH = '/home/somayeh/miniconda3/bin/spades.py'
#PATH_PREFIX = '/media/Data/PostDoc/Dalhousie/Work/Test1/'
PATH_PREFIX = '/media/Data/PostDoc/Dalhousie/Work/Test2/'
#PATH_PREFIX = '/media/Data/PostDoc/Dalhousie/Work/'

multi_processor = True
core_num = 4
coverage_thr = 30
task = [6]

#output_dir = 'Experiments/2_2_2/'
#output_dir = 'Experiments/real_sample1/'
output_dir = 'Experiments/CAMI_H_1/'
#output_dir = 'Experiments/metacherchant/CAMI_H_1/'
#output_dir = 'Experiments/northwood/'
#output_dir = 'Experiments/spacegraphcats/2_2_2/'

##amr_files = PATH_PREFIX + output_dir + 'AMR/sequences/TEM-181.fasta'
amr_files = PATH_PREFIX + output_dir +'AMR_info/sequences/'
artificial_amr_insertion = False
find_amr_genes = False
amr_identity_threshold =  95

ref_genomes_available = True
#ref_ng_seqs_file = PATH_PREFIX + output_dir + 'AMR_info/AMR_ref_neighborhood.fasta'
ref_ng_annotations_file = PATH_PREFIX + output_dir + 'AMR_info/ref_neighborhood_annotations.csv'
#ref_genome_files = PATH_PREFIX + output_dir + 'metagenome_data/'
#ref_genome_files = [PATH_PREFIX + 'SE_FDAARGOS_768.fasta', PATH_PREFIX + 'enterococcus_zy2.fasta']
#ref_genome_files = PATH_PREFIX + output_dir + 'M2_S001__insert_180_gsa_anonymous.fasta'
#ref_genome_files = PATH_PREFIX + output_dir + 'M2_S002__insert_180_gsa_anonymous.fasta'
ref_genome_files = PATH_PREFIX + output_dir + 'H_S001__insert_180_gsa_anonymous.fasta'


#insertion of AMR gene in ref_genome
number_of_insertions = [2, 1]
number_of_copies = [[3, 1],[1]]
insertion_type = Insertion_type.random
insertion_locations = []

genome_amr_files = [PATH_PREFIX + output_dir + "SE_FDAARGOS_768_TRU-1_37589.fasta",
                    PATH_PREFIX + output_dir + "enter_zy2_TRU1_19098.fasta"]

#simulating reads
read_length =  150
metagenome_file = PATH_PREFIX + output_dir +'metagenome.fasta'
#reads =[PATH_PREFIX + output_dir + 'sub50_trimmed_ERR1713331_1.fastq',\
# 		PATH_PREFIX + output_dir + 'sub50_trimmed_ERR1713331_2.fastq']
#reads = PATH_PREFIX + output_dir +'H_S001__insert_180_reads_anonymous.fq.gz'
reads =[PATH_PREFIX + output_dir + 'NW016_1_trimmed._R1.fastq',\
 		PATH_PREFIX + output_dir + 'NW016_1_trimmed._R2.fastq']

spades_thread_num = 16

assembler = Assembler_name.metaspades
#Setting for assembler
if assembler == Assembler_name.metaspades:
	assembler_output_dir = 'spade_output'
	gfa_file = PATH_PREFIX + output_dir + assembler_output_dir + '/assembly_graph_with_scaffolds.gfa'
	contig_file = PATH_PREFIX + output_dir + assembler_output_dir + '/contigs.fasta'
	max_kmer_size = 55
elif assembler == Assembler_name.megahit:
	assembler_output_dir = 'megahit_output'
	gfa_file = PATH_PREFIX + output_dir + assembler_output_dir + '/k59.gfa'
	#gfa_file = PATH_PREFIX + output_dir + assembler_output_dir + '/k119.gfa'
	#contig_file = PATH_PREFIX + output_dir + assembler_output_dir + '/final.contigs.fa'
	contig_file = PATH_PREFIX + output_dir + assembler_output_dir + '/intermediate_contigs/k59.contigs.fa'
	#max_kmer_size = 59
	max_kmer_size = 119
elif assembler == Assembler_name.bcalm:
	assembler_output_dir = 'bcalm_output'
	gfa_file = PATH_PREFIX + output_dir + assembler_output_dir + '/bcalm_graph_55.gfa'
	#doesn't produce any contig file; so just sue the one from meta-pades for evaluation
	contig_file = PATH_PREFIX + output_dir + 'spades_output/contigs.fasta'
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
#nodes in our path or already extracting path_seq_len_percent_threshod percent of
#neighborhood sequence, we stop or only traverse the path from the longest neighbor of the current node
path_node_threshold =  10
path_seq_len_percent_threshod = 90

PROKKA_COMMAND_PREFIX = 'docker run -v `pwd`:/data staphb/prokka:latest '

use_RGI =  True
RGI_include_loose = False

#neighborhood_seq_files = PATH_PREFIX + output_dir + 'assembly_neighborhood_sequences_5000_2020-09-27_20-55.txt'
ng_seq_files = PATH_PREFIX + output_dir + 'sequences_info/sequences_info_'+str(seq_length)+'/sequences/'
#neighborhood_seq_files = PATH_PREFIX + output_dir + 'sequences_info/sequences/ng_sequences_TEM-7_1000_2020-11-20_12-19.txt'
ng_path_info_files = PATH_PREFIX + output_dir + 'sequences_info/sequences_info_'+str(seq_length)+'/paths_info/'
