
from full_pipeline import Pipeline_tasks
from full_pipeline import Insertion_type

BANDAGE_PATH = '/media/Data/tools/Bandage_Ubuntu_dynamic_v0_8_1/Bandage'
ART_PATH = '/media/Data/tools/art_src_MountRainier_Linux/art_illumina'
SPADES_PATH = '/home/somayeh/miniconda3/bin/spades.py'
PATH_PREFIX = '/media/Data/PostDoc/Dalhousie/Work/Test1/'

#task = [Pipeline_tasks.all.value]
task = [5, 6]

amr_file = PATH_PREFIX + 'TRU-1.fasta'

#ref_genome_number = 2
ref_genome_files = [PATH_PREFIX + 'SE_FDAARGOS_768.fasta', PATH_PREFIX + 'enterococcus_zy2.fasta']
#ref_genome_files = [PATH_PREFIX + 'SE_FDAARGOS_768.fasta']

#output_dir = 'Set33/'
output_dir = 'Experiments/Set30_13Sep2020/salmonella7_enterococcus1/'

#insertion of AMR gene in ref_genome
number_of_insertions = [2, 1]
number_of_copies = [[3, 1],[1]]
insertion_type = Insertion_type.random
insertion_locations = []

#genome_amr_files = [PATH_PREFIX + output_dir + "SE_FDAARGOS_768_TRU-1_25292.fasta",
#                    PATH_PREFIX + output_dir + "SE_FDAARGOS_768_TRU-1_49186.fasta"]
genome_amr_files = [PATH_PREFIX + output_dir + "SE_FDAARGOS_768_TRU-1_37589.fasta",
                    PATH_PREFIX + output_dir + "enter_zy2_TRU1_19098.fasta"]
#genome_amr_files = [PATH_PREFIX + output_dir + "SE_FDAARGOS_768_TRU-1_5008.fasta",
#                    PATH_PREFIX + output_dir + "enter_zy2_TRU1_24316.fasta"]
#genome_amr_files = [PATH_PREFIX + output_dir + "SE_FDAARGOS_768_TRU-1_17439.fasta",
#                    PATH_PREFIX + output_dir + "SE_FDAARGOS_768_TRU-1_51890.fasta"]

#simulating reads
read_length =  150

read1=PATH_PREFIX + output_dir + "metagenome_1.fq"
read2=PATH_PREFIX + output_dir + "metagenome_2.fq"
spades_thread_num = 16
#spades_output_dir = 'spades_output'
spades_output_dir = 'megahit_output'
spades_error_correction = False

#gfa_file = PATH_PREFIX + output_dir + spades_output_dir + '/assembly_graph_with_scaffolds.gfa'
#contig_file = PATH_PREFIX + output_dir + spades_output_dir + '/contigs.fasta'
gfa_file = PATH_PREFIX + output_dir + spades_output_dir + '/k99.gfa'
contig_file = PATH_PREFIX + output_dir + spades_output_dir + '/final.contigs.fa'
graph_distance = 3
seq_length = 1000

PROKKA_COMMAND_PREFIX = 'docker run -v `pwd`:/data staphb/prokka:latest '

use_RGI =  True
RGI_include_loose = False

neighborhood_seq_file = PATH_PREFIX + output_dir + 'assembly_neighborhood_sequences_5000_2020-09-27_20-55.txt'
