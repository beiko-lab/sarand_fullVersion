import sys
import os
import csv
import argparse
import datetime
import logging
import shutil
#import sarand
import params
from __init__ import __version__
from full_pipeline import create_arguments, modify_params, full_pipeline_main
from utils import initialize_logger
from find_amrs_in_sample import create_ref_arguments, find_ref_amrs_main
from amr_neighborhood_in_contigs import create_contig_arguments, find_contig_amrs_main

def full_pipeline_init(args, params):
    """
    """
    #create the output directory; if it exists, delete it and create a new one
    if not os.path.exists(args.output_dir):
        os.makedirs(args.output_dir)
    else:
        try:
            shutil.rmtree(args.output_dir)
        except OSError as e:
            logging.error("Error: %s - %s." % (e.filename, e.strerror))
        os.makedirs(args.output_dir)
    logging.info('Running full_pipeline ...')
    params = modify_params(params, args)
    full_pipeline_main(params)

def find_ref_amrs_init(args, params):
    """
    """
    #create the output directory; if it exists, delete it and create a new one
    if not os.path.exists(args.output_dir):
        os.makedirs(args.output_dir)
    else:
        try:
            shutil.rmtree(args.output_dir)
        except OSError as e:
            logging.error("Error: %s - %s." % (e.filename, e.strerror))
        os.makedirs(args.output_dir)
    logging.info("Running find_amrs_in_sample ...")
    find_ref_amrs_main(args)

def find_contig_amrs_init(args, params):
    """
    """
    logging.info("Running amr_neighborhood_in_contigs ...")
    find_contig_amrs_main(args)

def main():
    parser = argparse.ArgumentParser(description="extract the neighborhood of the "
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
                                       usage="sarand_main.py full_pipeline <options>",
                                       help='')
    full_parser = create_arguments(params, full_parser)
    full_parser.set_defaults(func = full_pipeline_init)
    # add subparser for find_amrs_in_sample.py
    ref_parser = subparser.add_parser('find_ref_amrs', description="To find all AMRs available "
                                        "in a metagenome sample, extract their neighborhood "
                                        "sequences and annotate them",
                                       usage="sarand_main.py find_ref_amrs <options>",
                                       help='')
    ref_parser = create_ref_arguments(params, ref_parser)
    ref_parser.set_defaults(func = find_ref_amrs_init)
    # add subparser for amr_neighborhood_in_contigs.py
    contig_parser = subparser.add_parser('find_contig_amrs', description="To find the neighborhood "
                                            "of AMRs in a contig file, compare them with "
    			                            "that of the ref genomes and calculate the "
                                            "sentivity and precision",
                                       usage="sarand_main.py find_contig_amrs <options>",
                                       help='')
    contig_parser = create_contig_arguments(params, contig_parser)
    contig_parser.set_defaults(func = find_contig_amrs_init)

    args = parser.parse_args()
    #logging file
    log_name = 'logger_'+datetime.datetime.now().strftime('%Y-%m-%d_%H-%M')+'.log'
    initialize_logger(args.main_dir, log_name)
    logging.info(str(params.__dict__))
    args.func(args,params)


if __name__ == '__main__':
    main()
