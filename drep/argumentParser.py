#!/usr/bin/env python3

'''
dRep- parse command-line arguemnts
'''

__author__ = "Matt Olm"
__license__ = "MIT"
__email__ = "mattolm@gmail.com"
__status__ = "Development"

import argparse
import os
import sys

import drep
from drep.controller import Controller

def version():
    versionFile = open(os.path.join(drep.__path__[0], 'VERSION'))
    return versionFile.read().strip()

VERSION = version()

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

def printHelp():
    print('')
    print('                ...::: dRep v' + VERSION + ' :::...''')
    print('''\

  Choose one of the operations below for more detailed help.
  Example: dRep dereplicate_wf -h

  Workflows:
    dereplicate_wf  -> Combine several of the operations below to de-replicate a genome list
    compare_wf      -> Simply compare a list of genomes

  Single operations:
    filter          -> Filter a genome list based on size, completeness, and/or contamination
    cluster         -> Compare and cluster a genome list based on MASH and ANIn/gANI
    adjust          -> Adjust genome clusters
    choose          -> Choose the best genome from each genome cluster
    evaluate        -> Evaluate genome de-replication
    bonus           -> Other random operations (currently just determine taxonomy)
    analyze         -> Make figures related to the above operations; test alternative clustering
    ''')
def parse_args(args):
    parser = argparse.ArgumentParser(formatter_class=SmartFormatter)
    subparsers = parser.add_subparsers(help='Desired operation',dest='operation')

    #
    # Make a parent parser for all of the subparsers
    #
    parent_parser = argparse.ArgumentParser(add_help=False)
    parent_parser.add_argument("work_directory",help="R|Directory where data and output\
    \n*** USE THE SAME WORK DIRECTORY FOR ALL DREP OPERATIONS ***")

    Bflags = parent_parser.add_argument_group('SYSTEM PARAMETERS')
    Bflags.add_argument('-p','--processors',help='threads',default=6,type=int)
    Bflags.add_argument('-d','--dry',help='dry run- dont do anything',default=False,
                        action= "store_true")
    Bflags.add_argument('-o','--overwrite',help='overwrite existing data in work folder',
                        default=False, action= "store_true")
    Bflags.add_argument("-h", "--help", action="help", help="show this help message and exit")

    #
    # Make a parent parser for filter operation
    #
    filter_parent = argparse.ArgumentParser(add_help=False)
    fiflags = filter_parent.add_argument_group('FILTERING OPTIONS')
    fiflags.add_argument("-l","--length", help= "Minimum genome length",default=50000,
                            type = float)
    fiflags.add_argument("-comp","--completeness", help="Minumum genome completeness",
                            default = 75, type = float)
    fiflags.add_argument("-con","--contamination", help="Maximum genome contamination",
                            default = 25, type = float)
    fiflags.add_argument("-str","--strain_htr", help="Maximum strain heterogeneity",
                            default = 25, type = float)
    fiflags.add_argument("--skipCheckM", help="Don't run checkM- will ignore con and "\
                            + "comp settings", action='store_true')

    #
    # Make a parent parser for the cluster operation
    #
    cluster_parent = argparse.ArgumentParser(add_help=False)
    # Comparison Parameters
    Compflags = cluster_parent.add_argument_group('CLUSTERING PARAMETERS')
    Compflags.add_argument("-ms","--MASH_sketch",help="MASH sketch size", default=1000)
    Compflags.add_argument("-pa","--P_ani",help="ANI threshold to form primary (MASH) clusters",
                        default=0.9, type = float)
    Compflags.add_argument("--S_algorithm", help="Algorithm for secondary clustering comaprisons",
                        default='ANIn', choices={'ANIn','gANI'})
    Compflags.add_argument("-sa", "--S_ani", help="ANI threshold to form secondary clusters",
                        default=0.99, type = float)
    Compflags.add_argument("-nc", "--cov_thresh", help="Minmum level of overlap between\
        genomes when doing secondary comparisons", default=0.1)
    Compflags.add_argument("-cm", "--coverage_method", help="R|Method to calculate coverage of an alignment\n" \
        + "(for ANIn only; gANI can only do larger method)\n"
        + "total   = 2*(aligned length) / (sum of total genome lengths)\n" \
        + "larger  = max((aligned length / genome 1), (aligned_length / genome2)",
                        choices=['total', 'larger'], default='total')
    Compflags.add_argument("-n_PRESET", help= "R|Presets to pass to nucmer\n" \
        + "tight   = only align highly conserved regions\n" \
        + "normal  = default ANIn parameters", choices=['normal','tight'],default='normal')
    Compflags.add_argument("--clusterAlg", help="Algorithm used to cluster genomes (passed\
                        to scipy.cluster.hierarchy.linkage",default='average')
    Compflags.add_argument("--SkipMash",help="Skip MASH clustering,\
                        just do secondary clustering on all genomes",action='store_true')
    Compflags.add_argument("--SkipSecondary", help="Skip secondary clustering, just perform MASH\
                        clustering", action='store_true')

    #
    # Make a parent parser for scoring
    #
    scoring_parent = argparse.ArgumentParser(add_help=False)
    Sflags = scoring_parent.add_argument_group('SCORING CRITERIA\n'+
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

    #
    # Make a parent parser for taxonomy
    #
    tax_parent = argparse.ArgumentParser(add_help=False)
    Tflags = tax_parent.add_argument_group('TAXONOMY')
    Tflags.add_argument("--run_tax",help='generate taxonomy information (Tdb)', \
                    action = "store_true")
    Tflags.add_argument("--cent_index",help='path to centrifuge index (for example, ' + \
                    "/home/mattolm/download/centrifuge/indices/b+h+v")

    #
    # Make a parent parser for evaluate
    #
    evaluate_parent = argparse.ArgumentParser(add_help=False)
    Fflags = evaluate_parent.add_argument_group('WARNINGS')
    Fflags.add_argument("--warn_dist", default= 0.25, help= "How far from the threshold " +\
                        " to throw cluster warnings")
    Fflags.add_argument("--warn_sim", default= 0.98, help= "Similarity threshold for " +\
                        " warnings between dereplicated genomes")
    Fflags.add_argument("--warn_aln", default= 0.25, help= "Minimum aligned fraction for " +\
                        " warnings between dereplicated genomes (ANIn)")

    #
    # Make a parent parser for adjust
    #
    adjust_parent = argparse.ArgumentParser(add_help=False)
    Nflags = adjust_parent.add_argument_group('RE-CLUSTER PRIMARY CLUSETERS')
    Nflags.add_argument("-c","--cluster",help='primary cluster to be adjusted')
    Nflags.add_argument("-t",'--threshold',help= 'clustering threshold to apply',default=\
                        .99)
    Nflags.add_argument("-m",'--clustering_method',help= 'Clustering method to apply',\
                        choices = {'ANIn','gANI'}, default='ANIn')
    Nflags.add_argument("-mc",'--minimum_coverage',help= 'Minimum coverage for ANIn',\
                        default=0.1)
    Nflags.add_argument("-a","--clusterAlg", help="Algorithm used to cluster genomes (passed\
                        to scipy.cluster.hierarchy.linkage)",default='average',choices=\
                        {'single','complete','average','weighted'})

    '''
    ####### Arguments for filter operation ######
    '''

    filter_parser = subparsers.add_parser("filter",formatter_class=SmartFormatter,\
                    parents = [parent_parser, filter_parent], add_help=False)

    # I/O Parameters
    Iflags = filter_parser.add_argument_group('I/O PARAMETERS')
    Iflags.add_argument('-g','--genomes',nargs='*',help='genomes to filter in .fasta format.\
                        Not necessary if Bdb or Wdb already exist')
    Iflags.add_argument('--Chdb',help='checkM run already completed. Must be in \
                        --tab_table format.')
    Iflags.add_argument("--checkM_method", help="Either lineage_wf (more accurate) "\
                            + "or taxonomy_wf (faster)", choices={'taxonomy_wf','lineage_wf'},\
                            default = 'lineage_wf')

    '''
    ####### Arguments for clustering operation ######
    '''

    cluster_parser = subparsers.add_parser("cluster",formatter_class=SmartFormatter,\
                    parents = [parent_parser, cluster_parent], add_help=False)

    # I/O Parameters
    Iflags = cluster_parser.add_argument_group('I/O PARAMETERS')
    Iflags.add_argument('-g','--genomes',nargs='*',help='genomes to cluster in .fasta format.\
                        Not necessary if already loaded sequences with the "filter" operation')

    '''
    ####### Arguments for choose operation ######
    '''

    choose_parser = subparsers.add_parser("choose",formatter_class=SmartFormatter,\
                    parents = [parent_parser, scoring_parent], add_help=False)

    # Other
    Oflags = choose_parser.add_argument_group('OTHER')
    Oflags.add_argument("--checkM_method", help="Either lineage_wf (more accurate) "\
                            + "or taxonomy_wf (faster)", choices={'taxonomy_wf','lineage_wf'},\
                            default = 'lineage_wf')

    '''
    ####### Arguments for analyze operation ######
    '''

    analyze_parser = subparsers.add_parser("analyze",formatter_class=SmartFormatter,\
                    parents = [parent_parser, adjust_parent], add_help=False)

    # Plotting
    Caflags = analyze_parser.add_argument_group('PLOTTING')
    Caflags.add_argument("-pl", "--plots", help= "R|Plots. "
                        + "Input 'all' or 'a' to plot all\n"
                        + "1) Primary clustering dendrogram\n"
                        + "2) Secondary clustering dendrograms\n"
                        + "3) Secondary clusters heatmaps\n"
                        + "4) Comparison scatterplots\n"
                        + "5) Cluster scorring plot\n"
                        + "6) Winning genomes\n",
                        nargs='*')

    '''
    ####### Arguments for adjust operation ######
    '''

    adjust_parser = subparsers.add_parser("adjust",formatter_class=SmartFormatter,\
                    parents = [parent_parser, adjust_parent], add_help=False)

    # Remove clusters
    Rflags = adjust_parser.add_argument_group('CLUSTER REMOVAL')
    Rflags.add_argument("--rm_cluster",help='cluster(s) to be removed. Can be primary ' +\
                        "or secondary cluster(s). Will delete cluster from " +\
                        "Cdb, linkage (if primary cluster), Wdb, and ./dereplicated_genomes",
                        nargs='*')

    '''
    ####### Arguments for bonus operation ######
    '''

    bonus_parser = subparsers.add_parser("bonus",formatter_class=SmartFormatter,\
                    parents = [parent_parser, tax_parent], add_help=False)

    '''
    ####### Arguments for evaluate operation ######
    '''

    evaluate_parser = subparsers.add_parser("evaluate",formatter_class=SmartFormatter,\
                    parents = [parent_parser, evaluate_parent], add_help=False)

    # Evaluation
    Fflags = evaluate_parser.add_argument_group("EVALUATIONS")
    Fflags.add_argument("-e", "--evaluate", help= "R|Things to evaluate "
                        + "Input 'all' or 'a' to evaluate all\n"
                        + "1) Evaluate de-replicated genome similarity\n"
                        + "2) Throw warnings for clusters that were almost different\n"
                        + "3) Generate a database of information on winning genomes\n",
                        nargs='*')

    '''
    ####### Arguments for dereplicate_wf ######
    '''
    dereplicate_parser = subparsers.add_parser("dereplicate_wf",formatter_class=SmartFormatter,\
                        parents=[parent_parser, filter_parent, cluster_parent, scoring_parent,\
                        tax_parent, evaluate_parent], add_help=False, epilog=\
                        "Example: dRep dereplicate_wf output_dir/ -g /path/to/genomes/*.fasta")

    # I/O
    Iflags = dereplicate_parser.add_argument_group('I/O PARAMETERS')
    Iflags.add_argument('-g','--genomes',nargs='*',help='genomes to cluster in .fasta format')
    Iflags.add_argument("--checkM_method", help="Either lineage_wf (more accurate) "\
                            + "or taxonomy_wf (faster)", choices={'taxonomy_wf','lineage_wf'},\
                            default = 'lineage_wf')
    Iflags.add_argument('--Chdb',help='checkM run already completed. Must be in \
                        --tab_table format.')

    '''
    ####### Arguments for compare_wf ######
    '''
    dereplicate_parser = subparsers.add_parser("compare_wf",formatter_class=SmartFormatter,\
                        parents=[parent_parser, cluster_parent, tax_parent, evaluate_parent],\
                         add_help=False, epilog=\
                        "Example: dRep compare_wf output_dir/ -g /path/to/genomes/*.fasta")

    # I/O
    Iflags = dereplicate_parser.add_argument_group('I/O PARAMETERS')
    Iflags.add_argument('-g','--genomes',nargs='*',help='genomes to cluster in .fasta format')

    '''
    ####### PARSE THE ARGUMENTS ######
    '''

    # Handle the situation where the user wants the raw help
    #args = None
    #if (len(sys.argv) == 1 or sys.argv[1] == '-h' or sys.argv == '--help'):
    if (len(args) == 0 or args[0] == '-h' or args[0] == '--help'):
        printHelp()
        sys.exit(0)
    else:
        return parser.parse_args(args)
