#!/usr/bin/env python3

import argparse
import logging
import os

# !!! This is just for testing purposes, obviously
import sys
sys.path.append('/home/mattolm/Programs/drep/')
import drep_modules as dm
import drep_modules.WorkDirectory
import drep_modules.d_cluster

def drep_wrapper(args):
    """
    This is the meat of the program
    """
    
    # Make the WorkDirectory if it doesn't yet exist, and make a logger in it.
    args.work_directory = str(os.path.abspath(args.work_directory))
    if not os.path.exists(args.work_directory):
        os.makedirs(args.work_directory)
        
    log_dir = args.work_directory + '/log/'
    if not os.path.exists(log_dir):
        os.makedirs(log_dir)
        
    logging.basicConfig(filename=log_dir + 'logger.log',level=logging.DEBUG,\
        format='%(asctime)s %(message)s')
    logging.info("***Logger made at {0}***".format(args.work_directory + '/logger.log'))
    logging.info("Arguments: {0}".format(args))
    
    # Call the intended operation
    if args.operation == "cluster":
        cluster_operation(**vars(args))

def cluster_operation(**kwargs):
    
    logging.info("Starting the clustering operation")
    drep_modules.d_cluster.d_cluster_wrapper(kwargs['work_directory'],**kwargs)
    logging.info("Finished the clustering operation")

"""
########################################
#    Argument Parsing                  #
########################################
"""

class SmartFormatter(argparse.ArgumentDefaultsHelpFormatter):
    def _split_lines(self, text, width):
        if text.startswith('R|'):
            return text[2:].splitlines()  
        # this is the RawTextHelpFormatter._split_lines
        return argparse.HelpFormatter._split_lines(self, text, width)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=SmartFormatter)
    subparsers = parser.add_subparsers(help='Desired operation',dest='operation')
    
    ####### Arguments for clustering operation ######
    
    cluster_parser = subparsers.add_parser("cluster",formatter_class=SmartFormatter)
    cluster_parser.add_argument("work_directory",help="R|Directory where data and output\
    \n*** USE THE SAME WORK DIRECTORY FOR ALL DREP OPERATIONS ***")
    
    # Clustering Parameters
    Cflags = cluster_parser.add_argument_group('CLUSTERING PARAMETERS')
    Cflags.add_argument("-ma","--MASH_ani",help="ANI threshold to form MASH clusters",
                        default=0.9)
    Cflags.add_argument("-na", "--ANIn_ANI", help="ANI threshold to form ANIn clusters",
                        default=0.99)
    Cflags.add_argument("--SkipMash",help="Skip MASH clustering,\
                        just do ANIn on all genomes",action='store_true')
    Cflags.add_argument("--SkipANIn", help="Skip ANIn calculations, just perform MASH\
                        clustering", action='store_true')
        
    # Comparison Parameters
    Compflags = cluster_parser.add_argument_group('COMPARISON PARAMETERS')
    Compflags.add_argument("-ms","--MASH_sketch",help="MASH sketch size", default=1000)
    Compflags.add_argument("-nc", "--ANIn_cov", help="Fraction of genomes compared\
        required to form ANIn clusters", default=0.5)
    Compflags.add_argument("-n_c", help="nucmer minimum cluster length", default=60)
    Compflags.add_argument("-n_maxgap", help="nucmer maximum gap between two adjacent matches \
                        in a cluster", default=90)
    Compflags.add_argument("-n_noextend", help="toggle the nucmer cluster extension step", 
                        default=False, action="store_true")
    Compflags.add_argument("-n_method", help="nucmer method", default="mum")
    Compflags.add_argument("-n_PRESET", help=
        """R|*This will overwrite the above nucmer settings*
        tight   = only align highly conserved regions
        normal  = default ANIn parameters""")
        
    # I/O Parameters
    Iflags = cluster_parser.add_argument_group('I/O PARAMETERS')
    Iflags.add_argument('-g','--genomes',nargs='*',help='genomes to cluster in .fasta format.\
                        Not necessary if already loaded sequences with the "filter" operation')
        
    # Biotite Parameters
    Bflags = cluster_parser.add_argument_group('SYSTEM PARAMETERS')
    Bflags.add_argument('-p','--processors',help='threads',default=6)
    Bflags.add_argument('-d','--dry',help='dry run- dont do anything',default=False,
                        action= "store_true")
    Bflags.add_argument('-o','--overwrite',help='overwrite existing data in work folder',
                        default=False, action= "store_true")
    
    ####### Arguments for filtering operation ######
    
    filter_parser = subparsers.add_parser("filter")
    filter_parser.add_argument("-l",'--length',help="Minimum length")
    
    args = parser.parse_args()
    drep_wrapper(args)