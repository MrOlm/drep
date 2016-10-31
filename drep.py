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
import drep_modules.d_analyze
import drep_modules.d_filter
import drep_modules.d_choose

def drep_wrapper(args):
    """
    This is the meat of the program
    """
    # Call the intended operation
    if args.operation == "filter":
        filter_operation(**vars(args))
    if args.operation == "cluster":
        cluster_operation(**vars(args))
    if args.operation == "analyze":
        analyze_operation(**vars(args))
    if args.operation == "choose":
        choose_operation(**vars(args)) 

def cluster_operation(**kwargs):
    if (kwargs['MASH_ani'] > 1) or (kwargs['ANIn_ani'] > 1):
        print("Can't assign a MASH or ANIn value over 1")
        sys.exit()
    
    makeload_logger(kwargs['work_directory'])
    
    kwargs['M_Lcutoff'] = 1 - kwargs['MASH_ani']
    kwargs['M_clusterAlg'] = 'hierarchical'
    kwargs['M_Lmethod'] = kwargs['clusterAlg']
    
    kwargs['N_Lcutoff'] = 1 - kwargs['ANIn_ani']
    kwargs['N_clusterAlg'] = 'hierarchical'
    kwargs['N_Lmethod'] = kwargs['clusterAlg']
    
    logging.info("Starting the clustering operation")
    drep_modules.d_cluster.d_cluster_wrapper(kwargs['work_directory'],**kwargs)
    logging.info("!!! Finished the clustering operation !!!")

def analyze_operation(**kwargs):
    makeload_logger(kwargs['work_directory'])
    
    logging.info("Starting the analyze operation")
    drep_modules.d_analyze.d_analyze_wrapper(kwargs['work_directory'],**kwargs)
    logging.info("!!! Finished the analyze operation !!!")
    
def filter_operation(**kwargs):
    makeload_logger(kwargs['work_directory'])
    
    logging.info("Starting the filter operation")
    drep_modules.d_filter.d_filter_wrapper(kwargs['work_directory'],**kwargs)
    logging.info("!!! Finished the filter operation !!!")
    
def choose_operation(**kwargs):
    makeload_logger(kwargs['work_directory'])
    
    logging.info("Starting the choose operation")
    drep_modules.d_choose.d_choose_wrapper(kwargs['work_directory'],**kwargs)
    logging.info("!!! Finished the choose operation !!!")
    
def makeload_logger(wd):
    wd = str(os.path.abspath(wd))
    if not os.path.exists(wd):
        os.makedirs(wd)
        
    log_dir = wd + '/log/'
    if not os.path.exists(log_dir):
        os.makedirs(log_dir)
        
    logging.basicConfig(filename=log_dir + 'logger.log',level=logging.DEBUG,\
        format='%(asctime)s %(message)s')
    logging.info("***Logger started up at {0}***".format(log_dir + 'logger.log'))
    

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
    
    '''
    ####### Arguments for filter operation ######
    '''
    
    filter_parser = subparsers.add_parser("filter",formatter_class=SmartFormatter)
    filter_parser.add_argument("work_directory",help="R|Directory where data and output\
    \n*** USE THE SAME WORK DIRECTORY FOR ALL DREP OPERATIONS ***")
    
    # Filtering options, which will either filter Bdb or Wdb
    fiflags = filter_parser.add_argument_group('Filtering Options')
    fiflags.add_argument("-l","--length", help= "Minimum genome length",default=500000,
                            type = float)
    fiflags.add_argument("-comp","--completeness", help="Minumum genome completeness",
                            default = 75, type = float)
    fiflags.add_argument("-con","--contamination", help="Maximum genome contamination",
                            default = 25, type = float)
    fiflags.add_argument("--skipCheckM", help="Don't run checkM- will ignore con and "\
                            + "comp settings", action='store_true')
    fiflags.add_argument("--checkM_method", help="Either lineage_wf (more accurate) "\
                            + "or taxonomy_wf (faster)", choices={'taxonomy_wf','lineage_wf'},\
                            default = 'lineage_wf')
                            
    # I/O Parameters
    Iflags = filter_parser.add_argument_group('I/O PARAMETERS')
    Iflags.add_argument('-g','--genomes',nargs='*',help='genomes to filter in .fasta format.\
                        Not necessary if Bdb or Wdb already exist')
    Iflags.add_argument('--Chdb',help='checkM run already completed. Must be in \
                        --tab_table format.')
        
    # Biotite Parameters
    Bflags = filter_parser.add_argument_group('SYSTEM PARAMETERS')
    Bflags.add_argument('-p','--processors',help='threads',default=6)
    Bflags.add_argument('-d','--dry',help='dry run- dont do anything',default=False,
                        action= "store_true")
    Bflags.add_argument('-o','--overwrite',help='overwrite existing data in work folder',
                        default=False, action= "store_true")
    
    '''
    ####### Arguments for clustering operation ######
    '''
    
    cluster_parser = subparsers.add_parser("cluster",formatter_class=SmartFormatter)
    cluster_parser.add_argument("work_directory",help="R|Directory where data and output\
    \n*** USE THE SAME WORK DIRECTORY FOR ALL DREP OPERATIONS ***")
    
    # Clustering Parameters
    Cflags = cluster_parser.add_argument_group('CLUSTERING PARAMETERS')
    Cflags.add_argument("-ma","--MASH_ani",help="ANI threshold to form MASH clusters",
                        default=0.9, type = float)
    Cflags.add_argument("-na", "--ANIn_ani", help="ANI threshold to form ANIn clusters",
                        default=0.99, type = float)
    Cflags.add_argument("--clusterAlg", help="Algorithm used to cluster genomes (passed\
                        to scipy.cluster.hierarchy.linkage",default='average')
    Cflags.add_argument("--SkipMash",help="Skip MASH clustering,\
                        just do ANIn on all genomes",action='store_true')
    Cflags.add_argument("--SkipANIn", help="Skip ANIn calculations, just perform MASH\
                        clustering", action='store_true')
        
    # Comparison Parameters
    Compflags = cluster_parser.add_argument_group('COMPARISON PARAMETERS')
    Compflags.add_argument("-ms","--MASH_sketch",help="MASH sketch size", default=1000)
    Compflags.add_argument("-nc", "--ANIn_cov", help="Minmum level of overlap between\
        genomes when doing ANIn comparisons", default=0.5)
    Compflags.add_argument("-n_PRESET", help=
        """R|tight   = only align highly conserved regions
        normal  = default ANIn parameters""", choices=['normal','tight'],default='normal')
    
    '''
    Compflags.add_argument("-n_c", help="nucmer minimum cluster length", default=60)
    Compflags.add_argument("-n_maxgap", help="nucmer maximum gap between two adjacent matches \
                        in a cluster", default=90)
    Compflags.add_argument("-n_noextend", help="toggle the nucmer cluster extension step", 
                        default=False, action="store_true")
    Compflags.add_argument("-n_method", help="nucmer method", default="mum")
    '''
        
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
    
    '''
    ####### Arguments for choose operation ######
    '''
    
    choose_parser = subparsers.add_parser("choose",formatter_class=SmartFormatter)
    choose_parser.add_argument("work_directory",help="R|Directory where data and output\
    \n*** USE THE SAME WORK DIRECTORY FOR ALL DREP OPERATIONS ***")
    
    """
    # Winner analysis
    Wflags = choose_parser.add_argument_group('WINNER CHRITERIA')

    Wflags.add_argument("-comW","--completeness_weight", help= "Weight of genmecompleteness
    """
    
    # Scoring
    Sflags = choose_parser.add_argument_group('SCORRING CHRITERIA\n'+
              "Based off of the formula: Completeness - Contamination + log(N50) + log(size)")
    
    Sflags.add_argument("-comW","--completeness_weight" , default = 1, type= float,
                        help='completeness weight')
    Sflags.add_argument("-conW","--contamination_weight", default = 5, type= float,
                        help='contamination weight')
    Sflags.add_argument("-N50W","--N50_weight", default = 0.5, type= float,
                        help='weight of log(genome N50)')
    Sflags.add_argument("-sizeW","--size_weight", default = 0, type= float,
                        help='weight of log(genome size)')
    Sflags.add_argument("-strW","--strain_heterogeneity_weight", default = 1, type= float,
                        help='strain heterogeneity weight')
    
    # Biotite Parameters
    Bflags = choose_parser.add_argument_group('SYSTEM PARAMETERS')
    Bflags.add_argument('-p','--processors',help='threads',default=6)
    Bflags.add_argument('-d','--dry',help='dry run- dont do anything',default=False,
                        action= "store_true")
    Bflags.add_argument('-o','--overwrite',help='overwrite existing data in work folder',
                        default=False, action= "store_true")
    
    
    '''
    ####### Arguments for analyze operation ######
    '''
    
    analyze_parser = subparsers.add_parser("analyze",formatter_class=SmartFormatter)
    analyze_parser.add_argument("work_directory",help="R|Directory where data and output\
    \n*** USE THE SAME WORK DIRECTORY FOR ALL DREP OPERATIONS ***")
    
    # Plotting
    Caflags = analyze_parser.add_argument_group('PLOTTING')
    Caflags.add_argument("-cviz", "--cluster_visualization", help= "R|Plots. "
                        + "Input 'all' or 'a' to plot all\n" 
                        + "1) MASH cluster-map\n"
                        + "2) MASH dendrogram\n"
                        + "3) ANIn cluster-maps\n"
                        + "4) ANIn dendrograms\n"
                        + "5) Comparison scatterplots\n"
                        + "6) Simple bin scorring\n"
                        + "7) Complex bin scorring\n",
                        nargs='*')
                        
    # Cluster testing
    Mflags = analyze_parser.add_argument_group('TEST CLUSTERING')
    Mflags.add_argument("-m",'--mash_cluster',help= 'MASH cluster to test clustering on')
    Mflags.add_argument("-t",'--threshold',help= 'Clustering threshold to apply')
    Mflags.add_argument("-c",'--clustering_method',help= 'Clustering method to apply',\
                        choices = {'ANIn','gANI'})
    Mflags.add_argument("--clusterAlg", help="Algorithm used to cluster genomes (passed\
                        to scipy.cluster.hierarchy.linkage)",default='average',choices=\
                        {'single','complete','average','weighted'})
                        
    '''
    ####### Arguments for adjust operation ######
    '''
    
    adjust_parser = subparsers.add_parser("adjust",formatter_class=SmartFormatter)
    adjust_parser.add_argument("work_directory",help="R|Directory where data and output\
    \n*** USE THE SAME WORK DIRECTORY FOR ALL DREP OPERATIONS ***")
    
    # ANIn clusters
    Nflags = adjust_parser.add_argument_group('MASH CLUSTER ADJUSTMENT')
    Nflags.add_argument("-c","--mash_cluster",help='mash cluster to be adjusted')
    
    # Biotite Parameters
    Bflags = adjust_parser.add_argument_group('SYSTEM PARAMETERS')
    Bflags.add_argument('-p','--processors',help='threads',default=6)
    Bflags.add_argument('-d','--dry',help='dry run- dont do anything',default=False,
                        action= "store_true")
    Bflags.add_argument('-o','--overwrite',help='overwrite existing data in work folder',
                        default=False, action= "store_true")
                        

    args = parser.parse_args()
    drep_wrapper(args)