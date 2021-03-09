#!/usr/bin/env python

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

  Matt Olm. MIT License. Banfield Lab, UC Berkeley. 2017 (last updated 2020)

  See https://drep.readthedocs.io/en/latest/index.html for documentation
  Choose one of the operations below for more detailed help. 
  
  Example: dRep dereplicate -h

  Commands:
    compare            -> Compare and cluster a set of genomes
    dereplicate        -> De-replicate a set of genomes
    check_dependencies -> Check which dependencies are properly installed
    ''')


def parse_args(args):
    parser = argparse.ArgumentParser(formatter_class=SmartFormatter)
    subparsers = parser.add_subparsers(help='Desired operation', dest='operation')

    # Make a parent parser for all of the subparsers
    parent_parser = argparse.ArgumentParser(add_help=False)
    parent_parser.add_argument("work_directory", help="R|Directory where data and output are stored\
    \n*** USE THE SAME WORK DIRECTORY FOR ALL DREP OPERATIONS ***")

    Bflags = parent_parser.add_argument_group('SYSTEM PARAMETERS')
    Bflags.add_argument('-p', '--processors', help='threads', default=6, type=int)
    Bflags.add_argument('-d', '--debug', help='make extra debugging output', default=False,
                        action="store_true")
    Bflags.add_argument("-h", "--help", action="help", help="show this help message and exit")

    # Make a parent parser for genome input
    genome_parser = argparse.ArgumentParser(add_help=False)
    Gflags = genome_parser.add_argument_group('GENOME INPUT')
    Gflags.add_argument('-g', '--genomes', nargs='*', help='genomes to filter in .fasta format.\
                                    Not necessary if Bdb or Wdb already exist. Can also input a text file with paths to genomes, which results in fewer OS issues than wildcard expansion')

    # Make a parent parser for filter operation
    filtering_parent = argparse.ArgumentParser(add_help=False)
    fiflags = filtering_parent.add_argument_group('GENOME FILTERING OPTIONS')
    fiflags.add_argument("-l", "--length", help="Minimum genome length", default=50000,
                         type=float)
    fiflags.add_argument("-comp", "--completeness", help="Minumum genome completeness",
                         default=75, type=float)
    fiflags.add_argument("-con", "--contamination", help="Maximum genome contamination",
                         default=25, type=float)

    quality_parent = argparse.ArgumentParser(add_help=False)
    Iflags = quality_parent.add_argument_group('GENOME QUALITY ASSESSMENT OPTIONS')
    Iflags.add_argument("--ignoreGenomeQuality", help="Don't run checkM or do any \
                quality filtering. NOT RECOMMENDED! This is useful for use with bacteriophages\
                or eukaryotes or things where checkM scoring does not work. Will only \
                choose genomes based on length and N50", action='store_true')
    Iflags.add_argument('--genomeInfo', help='location of .csv file containing quality \
                    information on the genomes. Must contain: ["genome"(basename of .fasta file \
                    of that genome), "completeness"(0-100 value for completeness of the genome), \
                    "contamination"(0-100 value of the contamination of the genome)]')
    Iflags.add_argument("--checkM_method", help="Either lineage_wf (more accurate) \
                        or taxonomy_wf (faster)", choices={'taxonomy_wf', 'lineage_wf'}, \
                        default='lineage_wf')
    Iflags.add_argument("--set_recursion", help="Increases the python recursion limit. " \
                                                + "NOT RECOMMENDED unless checkM is crashing due to recursion issues. " \
                                                + "Recommended to set to 2000 if needed, but setting this could crash python", \
                        default='0')
    Iflags.add_argument("--checkm_group_size", help="The number of genomes passed to checkM at a time. Increasing this increases RAM but makes checkM faster", \
                        default=2000, type=int)

    # Make a parent parser for the cluster operation
    cluster_parent = argparse.ArgumentParser(add_help=False)

    Clustflags = cluster_parent.add_argument_group('GENOME COMPARISON OPTIONS')
    Clustflags.add_argument("--S_algorithm", help="R|Algorithm for secondary clustering comaprisons:\n" \
                                                  + "fastANI = Kmer-based approach; very fast\n" \
                                                  + "ANImf   = (DEFAULT) Align whole genomes with nucmer; filter alignment; compare aligned regions\n" \
                                                  + "ANIn    = Align whole genomes with nucmer; compare aligned regions\n" \
                                                  + "gANI    = Identify and align ORFs; compare aligned ORFS\n" \
                                                  + "goANI   = Open source version of gANI; requires nsmimscan\n",
                            default='ANImf', choices={'ANIn', 'gANI', 'ANImf', 'goANI', 'fastANI'})
    Clustflags.add_argument("-ms", "--MASH_sketch", help="MASH sketch size", default=1000)
    Clustflags.add_argument("--SkipMash", help="Skip MASH clustering,\
                            just do secondary clustering on all genomes", action='store_true')
    Clustflags.add_argument("--SkipSecondary", help="Skip secondary clustering, just perform MASH\
                            clustering", action='store_true')
    Clustflags.add_argument("--n_PRESET", help="R|Presets to pass to nucmer\n" \
                                               + "tight   = only align highly conserved regions\n" \
                                               + "normal  = default ANIn parameters", choices=['normal', 'tight'],
                            default='normal')

    Compflags = cluster_parent.add_argument_group('GENOME CLUSTERING OPTIONS')
    Compflags.add_argument("-pa", "--P_ani", help="ANI threshold to form primary (MASH) clusters",
                           default=0.9, type=float)
    Compflags.add_argument("-sa", "--S_ani", help="ANI threshold to form secondary clusters",
                           default=0.99, type=float)
    Compflags.add_argument("-nc", "--cov_thresh", help="Minmum level of overlap between\
        genomes when doing secondary comparisons", default=0.1, type=float)
    Compflags.add_argument("-cm", "--coverage_method", help="R|Method to calculate coverage of an alignment\n" \
                                                            + "(for ANIn/ANImf only; gANI and fastANI can only do larger method)\n"
                                                            + "total   = 2*(aligned length) / (sum of total genome lengths)\n" \
                                                            + "larger  = max((aligned length / genome 1), (aligned_length / genome2))\n",
                           choices=['total', 'larger'], default='larger')
    Compflags.add_argument("--clusterAlg", help="Algorithm used to cluster genomes (passed\
                        to scipy.cluster.hierarchy.linkage", default='average',
                           choices={'single', 'complete', 'average', 'weighted', 'centroid', 'median', 'ward'})

    GRflags = cluster_parent.add_argument_group('GREEDY CLUSTERING OPTIONS\n'
                                                'These decrease RAM use and runtime at the expense of a minor loss in '
                                                'accuracy.\nRecommended when clustering 5000+ genomes')
    GRflags.add_argument("--multiround_primary_clustering",
                         help='Cluster each primary clunk separately and '
                              'merge at the end with single linkage. Decreases '
                              'RAM usage and increases speed, and the cost of a minor '
                              'loss in precision and the inability to plot '
                              'primary_clustering_dendrograms. Especially helpful '
                              'when clustering 5000+ genomes. Will be done with '
                              'single linkage clustering',
                         action='store_true')
    GRflags.add_argument("--primary_chunksize",
                         help="Impacts multiround_primary_clustering. If you have more than this many genomes, "
                              "process them in chunks of this size.",
                         default=5000, type=int)
    GRflags.add_argument("--greedy_secondary_clustering",
                         help="Use a heuristic to avoid pair-wise comparisons when doing secondary clustering. Will "
                              "be done with single linkage clustering. Only works for fastANI S_algorithm option at "
                              "the moment",
                         action='store_true')
    GRflags.add_argument("--run_tertiary_clustering",
                            help="Run an additional round of clustering on the final genome set. This is especially "
                                 "useful when greedy clustering is performed and/or to handle cases where similar genomes "
                                 "end up in different primary clusters. Only works with dereplicate, not compare.",
                            action='store_true')

    # Make a parent parser for scoring
    scoring_parent = argparse.ArgumentParser(add_help=False)
    Sflags = scoring_parent.add_argument_group(
        "SCORING CRITERIA\nBased off of the formula: \nA*Completeness - B*Contamination + C*(Contamination * ("
        "strain_heterogeneity/100)) + D*log(N50) + E*log(size) + F*(centrality - S_ani)\n\nA = completeness_weight; B = "
        "contamination_weight; C = strain_heterogeneity_weight; D = N50_weight; E = size_weight; F = cent_weight")

    Sflags.add_argument("-comW", "--completeness_weight", default=1, type=float,
                        help='completeness weight')
    Sflags.add_argument("-conW", "--contamination_weight", default=5, type=float,
                        help='contamination weight')
    Sflags.add_argument("-strW", "--strain_heterogeneity_weight", default=1, type=float,
                        help='strain heterogeneity weight')
    Sflags.add_argument("-N50W", "--N50_weight", default=0.5, type=float,
                        help='weight of log(genome N50)')
    Sflags.add_argument("-sizeW", "--size_weight", default=0, type=float,
                        help='weight of log(genome size)')
    Sflags.add_argument("-centW", "--centrality_weight", default=1, type=float,
                        help='Weight of (centrality - S_ani)')
    Sflags.add_argument("-extraW", "--extra_weight_table", default=None, type=str,
                        help='Path to a tab-separated file with two-columns, no headers, listing genome and extra score to apply to that genome')

    # Make a parent parser for evaluate
    evaluate_parent = argparse.ArgumentParser(add_help=False)
    Fflags = evaluate_parent.add_argument_group('WARNINGS')
    Fflags.add_argument("--warn_dist", default=0.25, help="How far from the threshold " + \
                                                          " to throw cluster warnings")
    Fflags.add_argument("--warn_sim", default=0.98, help="Similarity threshold for " + \
                                                         " warnings between dereplicated genomes")
    Fflags.add_argument("--warn_aln", default=0.25, help="Minimum aligned fraction for " + \
                                                         " warnings between dereplicated genomes (ANIn)")

    dereplicate_parser = subparsers.add_parser("dereplicate", formatter_class=SmartFormatter, \
                                               parents=[parent_parser, genome_parser, filtering_parent, quality_parent,
                                                        cluster_parent, scoring_parent, \
                                                        evaluate_parent], add_help=False, epilog= \
                                                   "Example: dRep dereplicate output_dir/ -g /path/to/genomes/*.fasta")

    dereplicate_parser = subparsers.add_parser("compare", formatter_class=SmartFormatter, \
                                               parents=[parent_parser, genome_parser, cluster_parent, evaluate_parent], \
                                               add_help=False, epilog= \
                                                   "Example: dRep compare output_dir/ -g /path/to/genomes/*.fasta")

    dep_parser = subparsers.add_parser("check_dependencies", formatter_class=SmartFormatter)

    '''
    ####### PARSE THE ARGUMENTS ######
    '''

    # Handle the situation where the user wants the raw help
    if len(args) == 0 or args[0] == '-h' or args[0] == '--help':
        printHelp()
        sys.exit(0)
    else:
        return parser.parse_args(args)
