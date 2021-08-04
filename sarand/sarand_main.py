import sys
import os
import csv
import argparse
import datetime
import logging
import shutil
import yaml

from sarand import params, full_pipeline, utils, find_amrs_in_sample, amr_neighborhood_in_contigs
from sarand.__init__ import __version__
from sarand.full_pipeline import update_full_pipeline_params, full_pipeline_main
from sarand.utils import initialize_logger, validate_print_parameters_tools
from sarand.find_amrs_in_sample import find_ref_amrs_main, update_ref_params
from sarand.amr_neighborhood_in_contigs import find_contig_amrs_main, update_contig_params

def full_pipeline_init(args, params):
    """
    """
    # Rewrite parameters of params which are available in data and
    # are supposed to be set for this function
    params = update_full_pipeline_params(params, args)
    #logging file
    log_name = 'logger_'+datetime.datetime.now().strftime('%Y-%m-%d_%H-%M')+'.log'
    initialize_logger(params.main_dir, log_name)
    #logging.info(str(params.__dict__))

    #create the output directory;
    if not os.path.exists(params.output_dir):
        os.makedirs(params.output_dir)
    #if it exists, delete it and create a new one
    # else:
    #     try:
    #         shutil.rmtree(args.output_dir)
    #     except OSError as e:
    #         logging.error("Error: %s - %s." % (e.filename, e.strerror))
    #     os.makedirs(args.output_dir)
    logging.info('Running full_pipeline ...')
    params = validate_print_parameters_tools(params, "full_pipeline")
    full_pipeline_main(params)

def find_ref_amrs_init(args, params):
    """
    """
    # Rewrite parameters of params which are available in data and
    # are supposed to be set for this function
    params = update_ref_params(params, args)
    #logging file
    log_name = 'logger_'+datetime.datetime.now().strftime('%Y-%m-%d_%H-%M')+'.log'
    initialize_logger(params.main_dir, log_name)

    #create the output directory;
    if not os.path.exists(params.output_dir):
        os.makedirs(params.output_dir)
    #if it exists, delete it and create a new one
    # else:
    #     try:
    #         shutil.rmtree(args.output_dir)
    #     except OSError as e:
    #         logging.error("Error: %s - %s." % (e.filename, e.strerror))
    #     os.makedirs(args.output_dir)
    logging.info("Running find_amrs_in_sample ...")
    params = validate_print_parameters_tools(params, "find_ref_amrs")
    find_ref_amrs_main(params)

def find_contig_amrs_init(args, params):
    """
    """
    # Rewrite parameters of params which are available in data and
    # are supposed to be set for this function
    params = update_contig_params(params, args)
    #logging file
    log_name = 'logger_'+datetime.datetime.now().strftime('%Y-%m-%d_%H-%M')+'.log'
    initialize_logger(params.main_dir, log_name)
    #create the output directory;
    if not os.path.exists(params.output_dir):
        os.makedirs(params.output_dir)
    logging.info("Running amr_neighborhood_in_contigs ...")
    params = validate_print_parameters_tools(params, "find_contig_amrs")
    find_contig_amrs_main(params)

def main():
    parser = argparse.ArgumentParser(description="Extract the neighborhood of the "
                                        "target Antimicrobial Resistance (AMR) "
                                        "genes from the assembly graph.",
                                        prog='sarand',
                                        usage='sarand <tool> <options>')
    parser.add_argument('-v', '--version', action='version',
                        version=f"%(prog)s {__version__}")
    # add tool specific parsers
    subparser = parser.add_subparsers(title="Available tools under sarand",help='')
    # add subparser for full_pipeline.py
    full_parser = subparser.add_parser('full_pipeline', description="Complete pipeline "
                                        "to extract AMR neighborhood from assembly graph "
                                        "and annotate it",
                                       usage="sarand full_pipeline <options>",
                                       help='')
    #full_parser = create_arguments(params, full_parser)
    full_parser.add_argument('--config_file', '-C', type = str, default='',
		help = 'the config file to set parameters for full_pipeline()')
    full_parser.set_defaults(func = full_pipeline_init)

    # add subparser for find_amrs_in_sample.py
    ref_parser = subparser.add_parser('find_ref_amrs', description="To find all AMRs available "
                                        "in a metagenome sample, extract their neighborhood "
                                        "sequences and annotate them",
                                       usage="sarand find_ref_amrs <options>",
                                       help='')
    #ref_parser = create_ref_arguments(params, ref_parser)
    ref_parser.add_argument('--config_file', '-C', type = str, default='',
		help = 'the config file to set parameters for find_ref_amrs()')
    ref_parser.set_defaults(func = find_ref_amrs_init)

    # add subparser for amr_neighborhood_in_contigs.py
    contig_parser = subparser.add_parser('find_contig_amrs', description="To find the neighborhood "
                                            "of AMRs in a contig file, compare them with "
    			                            "that of the ref genomes and calculate the "
                                            "sensitivity and precision",
                                       usage="sarand find_contig_amrs <options>",
                                       help='')
    #contig_parser = create_contig_arguments(params, contig_parser)
    contig_parser.add_argument('--config_file', '-C', type = str, default='',
		help = 'the config file to set parameters for find_contig_amrs()')
    contig_parser.set_defaults(func = find_contig_amrs_init)

    args = parser.parse_args()
    #If no argument has been passed
    if not len(sys.argv) > 1:
        print("Please use -h option to access usage information!")
        sys.exit()
    # Check if the config file has correctly been set and exist???
    if args.config_file=='' or not os.path.isfile(args.config_file) or\
        not args.config_file.lower().endswith('yaml'):
        print(args.config_file+" doesn't exist or not in the right format! Please provide the correct path to the YAML config file!")
        sys.exit()
    # Read config file into a dictionery
    print("Reading the config file '"+args.config_file+"' ...")
    with open(args.config_file, 'r') as yamlfile:
        data = yaml.load(yamlfile, Loader=yaml.FullLoader)

    args.func(data, params)


if __name__ == '__main__':
    main()
