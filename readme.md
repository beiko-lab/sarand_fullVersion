# Documentation for AMR Context Extraction
This document provides an overview of how the proposed tool can be installed and used. It attempts to document all necessary details to set up the tool, and provides some run guides.

## Overview
This tool can be used to extract the neighborhood of the target Antimicrobial Resistance (AMR) genes from the assembly graph.
It can also be used to simulate sequence reads from some reference genomes (through ART), run MetaSPAdes to assemble the simulated reads and then reconstruct the neighborhood of the AMR genes.

## Installation
### Step I: install the dependencies
Our tool relies on several dependencies, including Python, Prokka, RGI, BLAST, Bandage, MetaSPAdes (in case the assembly graph has not already been generated) and ART (in case of simulating reads).
The most straight forward way to install this tool's dependencies is using bioconda.
#### Cloning the tool repository
`git clone https://github.com/beiko-lab/AMR_context`

Now, move to AMR_context directory.
#### Installing bioconda
Make sure [bioconda](https://bioconda.github.io/user/install.html) has been installed and the channels are set properly as follows.
```
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
```
#### Creating a new conda environment
It's recommend to set up a conda environment first, so packages and
versions don't get mixed up with system versions.
`conda create -n amr_context python=3.6.10`

#### Activating the conda environment
`conda activate amr_context`

#### Installing BLAST
The installation instructions are available for [Linux/Unix/MacOS](https://www.ncbi.nlm.nih.gov/books/NBK52640/) and [Windows](https://www.ncbi.nlm.nih.gov/books/NBK52637/). For more details, please refer to [Blast+ user manual](https://www.ncbi.nlm.nih.gov/books/NBK279690/).

#### Installing dependencies through conda
```
conda install pip
conda install rgi=5.1.1
conda install prokka
conda install art
conda install bandage
```
Note: In current implementation, Prokka has been run via docker. If you would like to install it from a different channel, you need to update PROKKA_COMMAND_PREFIX variable in params.py with an appropriate value which is probably an empty string.  

#### Installing python requirements
Note: Make sure that you are in the root directory of this tool (AMR_context).
`pip install -r requirements.txt`

### Step II: Testing
#### Updating params.py
Make sure to update the following parameters in code/params.py.
- BANDAGE_PATH (the path to access bandage executable file)
- ART_PATH (the path to access art_illumina directory)
- SPADES_PATH (the path to spades.py)
- PROKKA_COMMAND_PREFIX (it probably should either be an empty string or the command from docker if you installed prokka via docker)

To run the code, make sure you are in the created conda environment.
To activate it, run:

`conda activate amr_context`

and then run the code by:

`python code/full_pipeline.py`

Note: You don't need to install and set the parameters for ART and SPAdes if the assembly graph is provided as an input.

## Exploring the code
### Optional parameters to set
The list of all parameters that can be set in this tool have been provided in code/params.py, and includes the following parameters that can be set via command line as well:  
- -h, --help   
> show help message and exit
- --task TASK or [FIRST_TASK, LAST_TASK]  
> which task would you like to do? For the entire pipeline choose 0; otherwise either provide a number representing one of the following tasks or two numbers to denote the start and end tasks (and of course all tasks in the middle will be run too). Here is the list: metagenome_creation = 1 read_simulation = 2 assembly = 3 graph_neighborhood = 4 sequence_neighborhood = 5 neighborhood_annotation = 6, neighborhood_evaluation = 7
- --amr_file AMR_FILE, -A AMR_FILE  
> the path of the files containing the AMR genes sequence
- --ref_genome_files REF_GENOME_FILES [REF_GENOME_FILES ...]  
> the address of reference genomes; it can be a file, a list of files or a directory
- --main_dir MAIN_DIR, -M MAIN_DIR  
> the main dir to store all results
- --read_length READ_LENGTH  
> the length of simulated reads can be either 150 or 250  
- --spades_thread_num SPADES_THREAD_NUM  
> the number of threads used for MetaSPAdes  
- --spades_output_dir SPADES_OUTPUT_DIR  
> the output dir to store MetaSPAdes results  
- --graph_distance GRAPH_DISTANCE, -D GRAPH_DISTANCE  
> the maximum distance of neighborhood nodes to be extracted from the AMR gene  
- --seq_length SEQ_LENGTH, -L SEQ_LENGTH  
> the length of AMR gene's neighbourhood to be extracted  
- --ng_seq_files NEIGHBORHOOD_SEQ_FILE  
> the address of the files (directory) containing all extracted neighborhood sequences in assembly graph  
- --ng_path_info_files NG_PATH_INFO_FILES
> the address of the files (directory) containing all path information for extracted neighborhood sequences in assembly graph')
- --gfa_file GFA_FILE  
> the address of the file for assembly graph  
- --contig_file CONTIG_FILE  
> the address of the file containing contigs after assembly  
- --genome_amr_files GENOME_AMR_FILES [GENOME_AMR_FILES ...]  
> the address of the files containing genome after AMR insertion  
- --reads [READs]         
> the address of the files containing paired-end reads
- --spades_error_correction SPADES_ERROR_CORRECTION  
> Whether to turn on or off error correction in MetaSPAdes  
- --use_RGI USE_RGI     
> Whether to contribute RGI annotation in Prokka result  
- --RGI_include_loose RGI_INCLUDE_LOOSE  
> Whether to include loose cases in RGI result  
- --find_amr_genes <BOOLEAN>
> Whether to assume the AMR genes (in metagenome) are known or to look for them in assembly graph
- --amr_identity_threshold AMR_IDENTITY_THRESHOLD
> the threshold used for amr alignment: a hit is returned if identity/coverage >= threshold
- --path_node_threshold PATH_NODE_THRESHOLD
> the threshold used for recursive pre_path and post_path search as long as the length of the path is less that this threshold
- --path_seq_len_percent_threshold PATH_SEQ_LEN_PERCENT_THR
> the threshold used for recursive pre_seq and post_seq until we have this percentage of the required length after which we just extract from the longest neighbor
- --ref_genomes_available <BOOLEAN>
> Whether we have access to reference genome(s)
- --multi_processor <BOOLEAN>
> Whether to use multi processors for parallel programming
- --core_num CORE_NUM
> the number of cores used in case of parallel programming
- --coverage_thr COVERAGE_THRESHOLD
> coverage threshold to check if an annotated gene is truly AMR neighbor or just a false positive

### Python files
#### 1- full_pipeline.py
This is the core file to do all the steps available in our tool including concatenating ref genomes in a single file, simulating reads, assembling reads, extracting amr neighborhood, annotating amr neighborhood sequences and evaluation (in case that ref genomes are available).

To run, make sure that parameters are set in code/params.py:

`python code/full_pipeline.py`

#### 2- extract_neighborhood.py
This is the main file to extract the neighborhood of an AMR gene from an assembly graph.

To run:

```
python code/extract_neighborhood.py --amr/-A <AMR gene file path in FASTA format>
--gfa/-G <GFA assembly graph>
--length/-L <length of the linear sequence around AMR gene to be extracted (default = 1000)>
--main_dir <the output directory to store the results>
```
#### 3- find_amrs_in_sample.py
This code is used to find all AMRs available in a metagenome sample, extract their neighborhood sequences and annotate them.

To run:
```
python code/find_amrs_in_sample.py --db <metagenome file path>
    --seq <fasta file containing all AMR sequences>
```
Note: it reads 3 parameetrs from params.py:
- params.PROKKA_COMMAND_PREFIX
- params.use_RGI
- params.RGI_include_loose

#### 4- amr_neighborhood_in_contigs.py
This code is used to find the neighborhood of AMRs in a contig file, annotate them, compare them with that of the ref genomes and calculate the sentivity and precision.

To run:

`code/python amr_neighborhood_in_contigs.py`
NOTE: It reads required parameters from code/params.py and the most important parameters need to be set correctly there, are:
- params.seq_length
- params.contig_file
- params.amr_identity_threshold
- params.amr_files
- params.ref_ng_annotations_file
- params.main_dir

NOTE: The result are available in the following directory:
`params.main_dir+'contigs_output_'+str(params.seq_length)`

#### 5- annotation_visualization.py
This file is used to visualize sequences annotations.

To run:
```
python annotation_visualization.py --csvfile <annotation file path>
    --title <the image title> --output <the output image name>
```
