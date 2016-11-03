#!/usr/bin/env python3

import matplotlib
matplotlib.use('Agg')

import pandas as pd
import os
import seaborn as sns
from matplotlib import pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import scipy.cluster.hierarchy
import logging
import glob
import logging
import shutil
import math

import drep as dm
import drep
import drep.d_cluster as dClust
import drep.d_filter as dFilter

"""#############################################################################
                            MODULE ARCHITECTURE

****************************    Visualization   *******************************
Graphs - Scatterplots

*   plot_MASH_vs_ANIn_ani(Mdb, Ndb)
    - Plot MASH_ani vs. ANIn_ani (including correlation)

*   plot_MASH_vs_ANIn_cov(Mdb, Ndb)
    - Plot MASH_ani vs. ANIn_cov (including correlation)

*   plot_ANIn_vs_ANIn_cov(Mdb, Ndb)
    - Plot ANIn vs. ANIn_cov (including correlation)

*   plot_MASH_vs_len(Mdb, Ndb)
    - Plot MASH_ani vs. length_difference (including correlation)

*   plot_ANIn_vs_len(Ndb)
    - Plot ANIn vs. length_difference (including correlation)

Graphs - Cluster Visualization

*   plot_MASH_clusters(Mdb, linkage, threshold (optional))
    -  Make a dengrogram and a heatmap with clusters colored

*   plot_ANIn_clusters(Ndb, linkage, threshold (optional))
    - For each MASH cluster, make a dendrogram and heatmap with ANIn clusters colored

Graphs - Winner Visualization

*   plot_winner_scoring(Wdb)
    -  For each ANIn cluster, show the scoring

Graphs - Custom

*   plot_cluster_tightness(Ndb)
    - Come up with some way of visualizing the variation within ANIn clusters,
      versus the variation between clusters
      - Show the average and max tightness within and between all clusters

****************************    Clustering   *******************************

*   Idea here is that you try out a new clustering, and then decide if that's what you
    want to apply to adjust one of the clusters

################################################################################
"""
"""
WRAPPERS
"""

def d_analyze_wrapper(wd, **kwargs):

    # Load the workDirectory
    wd = drep.WorkDirectory.WorkDirectory(wd)

    if kwargs.get('plots') != None:
        logging.info("calling cluster_vis_wrapper")

        # If taxonomy info exists, add it
        if kwargs.get('include_taxonomy',True):
            Bdb = wd.get_db('Bdb')
            if 'taxonomy' in Bdb:
                genome2taxonomy = Bdb.set_index('genome')['taxonomy'].to_dict()
                kwargs['genome2taxonomy'] = genome2taxonomy

        cluster_vis_wrapper(wd, **kwargs)

    if kwargs.get('cluster') != None:
        logging.info("calling cluster_test_wrapper")
        cluster_test_wrapper(wd, **kwargs)

def cluster_vis_wrapper(wd, **kwargs):

    # Make the plot directory
    plot_dir = wd.location + '/figures/'
    if not os.path.exists(plot_dir):
        os.makedirs(plot_dir)

    # Determine what plots to make
    arg_plots = kwargs.get('plots')
    options = ['1','2','3','4','5','6']
    to_plot = parse_options(options, arg_plots)
    logging.info("making plots {0}".format(to_plot))
    print("making plots {0}".format(' '.join(sorted(to_plot))))

    # 1) Primary clustering dendrogram
    if '1' in to_plot:
        # Load the required data
        Mdb = wd.get_db('Mdb')
        Cdb = wd.get_db('Cdb')
        Pcluster = wd.get_primary_linkage()
        Plinkage = Pcluster['linkage']

        clust_args = wd.arguments['cluster']
        PL_thresh = clust_args.get('P_Lcutoff', False)

        # Make the plot
        print("Plotting primary dendrogram...")
        plot_MASH_dendrogram(Mdb, Cdb, Plinkage, threshold = PL_thresh,\
                        plot_dir = plot_dir)

    # 2) Secondary clustering dendrogram
    if '2' in to_plot:
        print("Plotting secondary dendrograms...")
        plot_secondary_dendrograms(wd, plot_dir, **kwargs)

    # 3) Secondary clusters heatmap
    if '3' in to_plot:
        # Load the required data
        print("Plotting secondary clusters heatmaps...")
        print("WRITE THIS PART")

    # 4) Comparison scatterplots
    if '4' in to_plot:
        # Load the required data
        Ndb = wd.get_db('Ndb')
        Mdb = wd.get_db('Mdb')

        # Make the plot
        print("Plotting Scatterplots...")
        plot_scatterplots(Mdb, Ndb, plot_dir = plot_dir)

    # 5) Simple bin scorring
    if '5' in to_plot:
        # Load the required data
        Sdb = wd.get_db('Sdb')
        Cdb = wd.get_db('Cdb')
        Wdb = wd.get_db('Wdb')

        # Make the plot
        print("Plotting simple scorring plot...")
        plot_winner_scoring_simple(Wdb, Sdb, Cdb, plot_dir = plot_dir)

    # 6) Complex bin scorring
    if '6' in to_plot:
        # Load the required data
        Sdb = wd.get_db('Sdb')
        Cdb = wd.get_db('Cdb')
        Wdb = wd.get_db('Wdb')
        Chdb = wd.get_db('Chdb')

        # Make the plot
        print("Plotting complex scorring plot...")
        plot_winner_scoring_complex(Wdb, Sdb, Cdb, Chdb, plot_dir = plot_dir, **kwargs)

def cluster_test_wrapper(wd, **kwargs):

    # Validate arguments
    cluster = kwargs.get('cluster')
    comp_method = kwargs.get('clustering_method','ANIn')
    assert comp_method in ['ANIn','gANI']
    clust_method = kwargs.get('clusterAlg')
    threshold = kwargs.pop('threshold',None)
    if threshold != None: threshold = 1- float(threshold)

    # Make a bdb listing the genomes to cluster
    Cdb = wd.get_db('Cdb')
    Bdb = wd.get_db('Bdb')
    genomes = Cdb['genome'][Cdb['primary_cluster'] == int(cluster)].tolist()
    bdb = Bdb[Bdb['genome'].isin(genomes)]


    # Make the comparison database
    Xdb = dClust.compare_genomes(bdb,comp_method,wd,**kwargs)

    # Make it symmetrical
    Xdb['av_ani'] = Xdb.apply(lambda row: dClust.average_ani (row,Xdb),axis=1)
    Xdb['dist'] = 1 - Xdb['av_ani']
    db = Xdb.pivot("reference","querry","dist")

    # Cluster it
    if threshold == None:
        threshold = float(0)
    cdb, linkage = dClust.cluster_hierarchical(db, linkage_method = clust_method, \
                            linkage_cutoff = threshold)

    # Make the plot directory
    plot_dir = wd.location + '/figures/cluster_tests/'
    if not os.path.exists(plot_dir):
        os.makedirs(plot_dir)

    # Make the plot
    names = list(db.columns)
    plot_single_dendrogram(linkage, names, threshold=threshold, \
        title="{0}_MASH_cluster_{1}".format(comp_method,cluster), loc = plot_dir, \
        **kwargs)


def parse_options(options, args):
    to_plot = []

    if args[0] in ['all','a']:
        to_plot += options

    else:
        for arg in args:
            if arg in options:
                to_plot.append(arg)
            else:
                for letter in arg:
                    if letter in options:
                        to_plot.append(letter)
                    else:
                        print("Can't interpret argument {0}- quitting".format(arg))
                        sys.exit()

    return to_plot

"""
HEAT MAPS
"""

def plot_Mdb_heatmap(Mdb):
    db = Mdb.pivot("genome1","genome2","similarity")
    g = sns.clustermap(db,method=METHOD)
    g.fig.suptitle("MASH ANI")
    plt.setp(g.ax_heatmap.yaxis.get_majorticklabels(), rotation=0)
    return g

def plot_ANIn_heatmap(Ndb):
    gs = []
    for Mcluster in Ndb['MASH_cluster'].unique():
        db = Ndb[Ndb['MASH_cluster'] == Mcluster]
        if len(db['reference'].unique()) == 1:
            continue
        d = db.pivot("reference","querry","ani")
        g = sns.clustermap(d,method=METHOD)
        g.fig.suptitle("MASH cluster {0} - ANIn".format(Mcluster))
        plt.setp(g.ax_heatmap.yaxis.get_majorticklabels(), rotation=0)
        gs.append(g)
    return gs

def plot_ANIn_cov_heatmap(Ndb):
    gs = []
    for Mcluster in Ndb['MASH_cluster'].unique():
        db = Ndb[Ndb['MASH_cluster'] == Mcluster].copy()
        if len(db['reference'].unique()) == 1:
            continue
        d = db.pivot("reference","querry","alignment_coverage")
        g = sns.clustermap(d,method=METHOD)
        g.fig.suptitle("MASH cluster {0} - Alignment Coverage".format(Mcluster))
        plt.setp(g.ax_heatmap.yaxis.get_majorticklabels(), rotation=0)
        gs.append(g)
    return gs

"""
SCATTER PLOTS
"""

def plot_scatterplots(Mdb, Ndb, plot_dir=False):
    save = False

    # Initialize a .pdf
    if plot_dir != False:
        pp = PdfPages(plot_dir + 'Clustering_scatterplots.pdf')
        save = True

    g = plot_MASH_vs_ANIn_ani(Mdb,Ndb,exclude_zero_MASH=False)
    if save: pp.savefig(g)
    plt.show()

    g = plot_MASH_vs_ANIn_cov(Mdb,Ndb)
    if save: pp.savefig(g)
    plt.show()

    g = plot_ANIn_vs_ANIn_cov(Ndb)
    if save: pp.savefig(g)
    plt.show()

    g = plot_MASH_vs_len(Mdb,Ndb)
    if save: pp.savefig(g)
    plt.show()

    g = plot_ANIn_vs_len(Mdb,Ndb)
    if save: pp.savefig(g)
    plt.show()

    if save:
        pp.close()
    plt.close('all')

def plot_MASH_vs_ANIn_ani(Mdb,Ndb,exclude_zero_MASH=True):
    plt.close('all')
    mdb = Mdb.copy()
    mdb.rename(columns={'genome1':'querry','genome2':'reference',
                        'similarity':'MASH_ANI'},inplace=True)
    if exclude_zero_MASH:
        mdb= mdb[mdb['MASH_ANI'] > 0]

    db = pd.merge(mdb,Ndb)
    db.rename(columns={'ani':'ANIn'},inplace=True)
    g = sns.jointplot(x='ANIn',y='MASH_ANI',data=db)
    return plt.gcf()

def plot_MASH_vs_ANIn_cov(Mdb,Ndb,exclude_zero_MASH=True):
    plt.close('all')
    mdb = Mdb.copy()
    mdb.rename(columns={'genome1':'querry','genome2':'reference',
                        'similarity':'MASH_ANI'},inplace=True)
    if exclude_zero_MASH:
        mdb= mdb[mdb['MASH_ANI'] > 0]

    db = pd.merge(mdb,Ndb)
    db.rename(columns={'alignment_coverage':'ANIn_alignment_coverage'},inplace=True)
    g = sns.jointplot(x='ANIn_alignment_coverage',y='MASH_ANI',data=db)
    return plt.gcf()

def plot_ANIn_vs_ANIn_cov(Ndb):
    plt.close('all')
    db = Ndb.copy()
    db.rename(columns={'alignment_coverage':'ANIn_alignment_coverage','ani':'ANIn'},inplace=True)
    g = sns.jointplot(y='ANIn',x='ANIn_alignment_coverage',data=db)
    return plt.gcf()

def plot_MASH_vs_len(Mdb,Ndb,exclude_zero_MASH=True):
    plt.close('all')
    mdb = Mdb.copy()
    mdb.rename(columns={'genome1':'querry','genome2':'reference',
                        'similarity':'MASH_ANI'},inplace=True)
    if exclude_zero_MASH:
        mdb= mdb[mdb['MASH_ANI'] > 0]

    db = pd.merge(mdb,Ndb)
    db['length_difference'] = abs(db['reference_length'] - db['querry_length'])
    g = sns.jointplot(x='MASH_ANI',y='length_difference',data=db)
    return plt.gcf()

def plot_ANIn_vs_len(Mdb,Ndb,exclude_zero_MASH=True):
    plt.close('all')
    mdb = Mdb.copy()
    mdb.rename(columns={'genome1':'querry','genome2':'reference',
                        'similarity':'MASH_ANI'},inplace=True)
    if exclude_zero_MASH:
        mdb= mdb[mdb['MASH_ANI'] > 0]

    db = pd.merge(mdb,Ndb)
    db['length_difference'] = abs(db['reference_length'] - db['querry_length'])
    db.rename(columns={'ani':'ANIn'},inplace=True)
    g = sns.jointplot(x='ANIn',y='length_difference',data=db)
    return plt.gcf()

"""
CLUSETER PLOTS
"""

'''
def plot_MASH_clustermap(Mdb, Cdb, linkage, threshold = False, plot_dir = False, **kwargs):
    db = Mdb.pivot("genome1","genome2","similarity")
    names = list(db.columns)
    name2cluster = Cdb.set_index('genome')['MASH_cluster'].to_dict()
    colors = gen_color_list(names, name2cluster)
    name2color = gen_color_dictionary(names, name2cluster)

    # Make the clustermap
    g = sns.clustermap(db, row_linkage = linkage, col_linkage = linkage, \
                        row_colors = colors, col_colors = colors, vmin = 0.8, \
                        vmax = 1)
    g.fig.suptitle("MASH clustering")
    plt.setp(g.ax_heatmap.yaxis.get_majorticklabels(), rotation=0)

    # Adjust the figure size
    plt.subplots_adjust(bottom=0.3, right=0.7)
    fig = plt.gcf()
    fig.set_size_inches(x_fig_size(len(names)),x_fig_size(len(names)))

    #Save the figure
    if plot_dir != None:
        plt.savefig(plot_dir + 'Cviz_MASH_clustermap.pdf')
    plt.show()
    plt.close('all')
'''

def plot_MASH_dendrogram(Mdb, Cdb, linkage, threshold = False, plot_dir = False, **kwargs):
    db = Mdb.pivot("genome1","genome2","similarity")
    names = list(db.columns)
    name2cluster = Cdb.set_index('genome')['primary_cluster'].to_dict()
    name2color = gen_color_dictionary(names, name2cluster)

    # Make the dendrogram
    g = fancy_dendrogram(linkage,names,name2color,threshold=threshold)
    plt.title('MASH clustering')
    plt.xlabel('MASH Average Nucleotide Identity (ANI)')
    #plt.xlim([0,.4])

    # Adjust the figure size
    fig = plt.gcf()
    fig.set_size_inches(10,x_fig_size(len(names)))
    plt.subplots_adjust(left=0.3)

    # Adjust the labels
    plt.tick_params(axis='both', which='major', labelsize=8)
    axes = plt.gca()
    labels = axes.xaxis.get_majorticklocs()
    for i, label in enumerate(labels):
        labels[i] = (1 - float(label)) * 100
    axes.set_xticklabels(labels)


    # Save the figure
    if plot_dir != None:
        plt.savefig(plot_dir + 'Primary_clustering_dendrogram.pdf')
    plt.show()
    plt.close('all')

'''
def plot_ANIn_clustermaps(Ndb, Cdb, cluster2linkage, threshold = False, plot_dir = False, **kwargs):
    save = False

    # Initialize a .pdf
    if plot_dir != False:
        pp = PdfPages(plot_dir + 'Cviz_ANIn_clustermaps.pdf')
        save = True

    for cluster in sorted([int(x) for x in cluster2linkage.keys()]):
        cluster = str(cluster)
        logging.info("plotting ANIn cluster {0}".format(cluster))
        linkage = cluster2linkage[cluster]

        # Filter Ndb to just have the clusters of the linkage
        c_genomes = Cdb['genome'][Cdb['MASH_cluster'] == int(cluster)]
        db = Ndb[Ndb['reference'].isin(c_genomes)]
        db = db.pivot("reference","querry","ani")

        # Get the colors set up
        names = list(db.columns)
        name2cluster = Cdb.set_index('genome')['secondary_cluster'].to_dict()
        name2color = gen_color_dictionary(names, name2cluster)
        colors = [name2color[x] for x in names]

        # Make the clustermap
        g = sns.clustermap(db, row_linkage = linkage, col_linkage = linkage, \
                            row_colors = colors, col_colors = colors, vmin= 0.95,\
                            vmax= 1)

        g.fig.suptitle('ANI of MASH cluster {0}'.format(cluster))
        plt.setp(g.ax_heatmap.yaxis.get_majorticklabels(), rotation=0)

        # Adjust the figure size
        plt.subplots_adjust(bottom=0.3, right=0.7)
        fig = plt.gcf()
        fig.set_size_inches(10,x_fig_size(len(names)))
        fig.suptitle("MASH Cluster {0}".format(cluster), fontsize=20)

        # Save the clustermap
        if save == True:
            pp.savefig(fig)
        plt.show()
        plt.close(fig)

    pp.close()
    plt.close('all')
'''

def plot_ANIn_dendrograms(Ndb, Cdb, cluster2linkage, threshold = False, plot_dir = False, **kwargs):
    save = False

    # Initialize a .pdf
    if plot_dir != False:
        pp = PdfPages(plot_dir + 'Cviz_ANIn_dendrograms.pdf')
        save = True

    for cluster in sorted([int(x) for x in cluster2linkage.keys()]):
        cluster = str(cluster)
        logging.info("plotting ANIn cluster {0}".format(cluster))
        linkage = cluster2linkage[cluster]

        # Filter Ndb to just have the clusters of the linkage
        c_genomes = Cdb['genome'][Cdb['primary_cluster'] == int(cluster)]
        db = Ndb[Ndb['reference'].isin(c_genomes)]

        # Get the highest self-comparison
        self_thresh = 1 - db['ani'][db['reference'] == db['querry']].min()
        if self_thresh == float(0):  self_thresh = 1.0e-7 # Still draw if 0

        # Get the colors set up
        db = db.pivot("reference","querry","ani")
        names = list(db.columns)
        name2cluster = Cdb.set_index('genome')['secondary_cluster'].to_dict()
        colors = gen_color_list(names, name2cluster)
        name2color = gen_color_dictionary(names, name2cluster)

        # Make the dendrogram
        g = fancy_dendrogram(linkage,names,name2color,threshold=threshold,self_thresh =\
                            self_thresh)
        plt.title('MASH cluster {0}'.format(cluster))
        plt.xlabel('Average Nucleotide Identity (ANI)')
        if threshold != False:
            plt.xlim([0,2*threshold])

        # Adjust the figure size
        fig = plt.gcf()
        fig.set_size_inches(10,x_fig_size(len(names)))
        plt.subplots_adjust(left=0.5)

        # Adjust the labels
        plt.tick_params(axis='both', which='major', labelsize=8)
        axes = plt.gca()
        labels = axes.xaxis.get_majorticklocs()
        for i, label in enumerate(labels):
            labels[i] = (1 - float(label)) * 100
        axes.set_xticklabels(labels)

        # Save the clustermap
        if save == True:
            pp.savefig(fig)
        plt.show()
        plt.close(fig)

    pp.close()
    plt.close('all')

def plot_secondary_dendrograms(wd, plot_dir, **kwargs):
    save = False
    win = False

    # Initialize a .pdf
    if plot_dir != False:
        pp = PdfPages(plot_dir + 'Secondary_clustering_dendrograms.pdf')
        save = True

    # Load secondary databases if they exist (Gdb, Ndb)
    Sdb = {}
    if wd.hasDb('Ndb'):
        Sdb['ANIn'] = wd.get_db('Ndb')
    if wd.hasDb('Gdb'):
        Sdb['gANI'] = wd.get_db('Gdb')

    # Load winner database if it exists
    if wd.hasDb('Wdb'):
        Wdb = wd.get_db('Wdb')
        winners = Wdb['genome'].unique()
        win = True

    # For every cluster:
    Cdb = wd.get_db('Cdb')
    for cluster in sorted(Cdb['primary_cluster'].unique()):
        d = Cdb[Cdb['primary_cluster'] == cluster]

        # Skip if it's a singleton
        if len(d['genome'].unique()) == 1:
            continue

        # Load the linkage information
        linkI = wd.get_cluster("secondary_linkage_cluster_{0}".format(cluster))
        db = linkI['db']
        linkage = linkI['linkage']
        args = linkI['arguments']
        threshold = args['linkage_cutoff']
        alg = args['comparison_algorithm']
        clust_alg = args['linkage_method']

        # Get the colors set up
        names = list(db.columns)
        name2cluster = Cdb.set_index('genome')['secondary_cluster'].to_dict()

        # Handle the case where you deleted a secondary cluster
        for name in names:
            if name not in Cdb['genome'].tolist():
                name2cluster[name] = '0_0'

        name2color = gen_color_dictionary(names, name2cluster)

        # Get the highest self-comparison
        self_thresh = False
        if alg in Sdb:
            self_thresh = get_highest_self(Sdb[alg], names)

        # Make the dendrogram
        sns.set_style('whitegrid')
        g = fancy_dendrogram(linkage,names,name2color,threshold=threshold,self_thresh =\
                            self_thresh)

        # Add the title and subtitle
        title_string = 'Primary cluster {0}'.format(cluster)
        subtitle_string = "Comparison method: {0}    Clustering method: {1}".format(alg, clust_alg)
        plt.suptitle(title_string, y=1, fontsize=18)
        plt.title(subtitle_string, fontsize=10)


        plt.xlabel('Average Nucleotide Identity (ANI)')
        if threshold != False:
            plt.xlim([0,2*threshold])

        # Adjust the figure size
        fig = plt.gcf()
        fig.set_size_inches(10,x_fig_size(len(names)))
        plt.subplots_adjust(left=0.5)

        # Adjust the labels
        plt.tick_params(axis='both', which='major', labelsize=12)
        axes = plt.gca()
        labels = axes.xaxis.get_majorticklocs()
        for i, label in enumerate(labels):
            labels[i] = (1 - float(label)) * 100
        axes.set_xticklabels(labels)
        plt.gca().yaxis.grid(False)

        # Mark winning ones
        if win:
            ax = plt.gca()
            labels = [item.get_text() for item in ax.get_yticklabels()]
            for i, label in enumerate(labels):
                if label in winners: labels[i] = label + ' *'
            ax.set_yticklabels(labels)

        # Add taxonomy
        if kwargs.get('genome2taxonomy',False) != False:
            g2t = kwargs.get('genome2taxonomy')
            axes = plt.gca()
            labels = [item.get_text() for item in axes.get_yticklabels()]
            for i, label in enumerate(labels):
                labels[i] = "{0}\n{1}".format(label, g2t[label.replace(' *','')])
            axes.set_yticklabels(labels)

        # Save the file
        if save == True:
            pp.savefig(fig)
        plt.show()
        plt.close(fig)

    pp.close()
    plt.close('all')


def get_highest_self(db, genomes, min = 1.0e-7):
    d = db[db['reference'].isin(genomes)]
    self_thresh = 1 - d['ani'][d['reference'] == d['querry']].min()

    # Because 0s don't show up on the graph
    if self_thresh == float(0):
        self_thresh =  min

    return self_thresh

def plot_single_dendrogram(linkage, names, threshold=False, title=None, loc = None, **kwargs):
    '''
    names can be gotten like:
    db = db.pivot("reference","querry","ani")
    names = list(db.columns)
    '''

    # Make the dendrogram
    g = fancy_dendrogram(linkage,names,threshold=threshold)
    plt.title(title)
    plt.xlabel(kwargs.get('xlabel','ANI'))
    if threshold != False:
        plt.xlim([0,2*threshold])

    # Adjust the figure size
    fig = plt.gcf()
    fig.set_size_inches(10,x_fig_size(len(names)))
    plt.subplots_adjust(left=0.5)

    # Adjust the labels
    plt.tick_params(axis='both', which='major', labelsize=8)
    axes = plt.gca()
    labels = axes.xaxis.get_majorticklocs()
    for i, label in enumerate(labels):
        labels[i] = (1 - float(label)) * 100
    axes.set_xticklabels(labels)

    # Save the clustermap
    if loc != None:
        plt.savefig("{0}{1}.pdf".format(loc,title))
    plt.show()
    plt.close(fig)


"""
WINNER PLOTS
"""
def plot_winner_scoring_simple(Wdb, Sdb, Cdb, plot_dir= False, **kwargs):
    save = False

    # Initialize a .pdf
    if plot_dir != False:
        pp = PdfPages(plot_dir + 'Simple_cluster_scoring.pdf')
        save = True

    winners = list(Wdb['genome'].unique())
    for cluster in sorted(Cdb['secondary_cluster'].unique(), key=lambda x: comp_cluster(x)):
        d = Cdb[Cdb['secondary_cluster'] == cluster]
        d = d.merge(Sdb, how='left', on= 'genome')

        g = sns.barplot(data=d, x='genome', y='score')
        plt.title('Scoring of cluster {0}'.format(cluster))
        plt.ylabel('Score')
        plt.xticks(rotation=90)

        # Mark winning one
        labels = d['genome'].tolist()
        for i, label in enumerate(labels):
            if label in winners: labels[i] = label + ' *'
        axes = plt.gca()
        axes.set_xticklabels(labels)

        fig = plt.gcf()
        fig.set_size_inches(6,10)
        plt.subplots_adjust(bottom=0.5)

        if save == True:
            pp.savefig(fig)
        plt.show()
        plt.close(fig)

    if save:
        pp.close()
    plt.close('all')

def plot_winner_scoring_complex(Wdb, Sdb, Cdb, Chdb, plot_dir= False, **kwargs):
    save = False

    # Initialize a .pdf
    if plot_dir != False:
        pp = PdfPages(plot_dir + 'Complex_cluster_scoring.pdf')
        save = True

    # Figure out what you're going to show
    winners = list(Wdb['genome'].unique())
    Chdb['genome'] = Chdb['Bin Id']
    bars = kwargs.get('to_show',['score','N50 (scaffolds)','Completeness','Contamination','Genome size (bp)', 'Strain heterogeneity'])

    for cluster in sorted(Cdb['secondary_cluster'].unique(), key=lambda x: comp_cluster(x)):
        # Make a db for this cluster
        d = Cdb[Cdb['secondary_cluster'] == cluster]
        d = d.merge(Sdb, how='left', on= 'genome')
        d = d.merge(Chdb, how='left', on= 'genome')
        d = d[bars + ['genome']]

        # Make the normalize bar plot
        nd = normalize(d)
        db = pd.melt(nd, id_vars=['genome'], value_vars=bars)
        g = sns.barplot(data=db, y='genome', x='value', hue='variable')

        # Get a list of the un-normalized values
        x = pd.melt(d, id_vars=['genome'], value_vars=bars)
        vals = []
        for variable in x['variable'].unique():
            vals += [v for v in x['value'][x['variable'] == variable].tolist()]

        # Add un-normalized values to barplots
        i = 0
        for p in g.patches:
            g.annotate("{0:.1f}".format(vals[i]), (p.get_width(), p.get_y()+(p.get_height()/1.1) ), fontsize=8)
            i += 1

        plt.title('Scoring of cluster {0}'.format(cluster))
        plt.xlabel('Normalized Score')
        plt.legend(loc=(0,0))
        plt.tick_params(axis='both', which='major', labelsize=8)

        # Mark winning one
        labels = d['genome'].tolist()
        for i, label in enumerate(labels):
            if label in winners: labels[i] = label + ' *'
        axes = plt.gca()
        axes.set_yticklabels(labels)

        # Add taxonomy
        if kwargs.get('genome2taxonomy',False) != False:
            g2t = kwargs.get('genome2taxonomy')
            axes = plt.gca()
            labels = [item.get_text() for item in axes.get_yticklabels()]
            for i, label in enumerate(labels):
                labels[i] = "{0}\n{1}".format(label, g2t[label.replace(' *','')])
            axes.set_yticklabels(labels)

        fig = plt.gcf()
        fig.set_size_inches(12,x_fig_size(len(labels), factor=1))
        plt.subplots_adjust(left=0.5)

        if save == True:
            pp.savefig(fig)
        plt.show()
        plt.close(fig)

    if save:
        pp.close()
    plt.close('all')

"""
OTHER
"""

def x_fig_size(points, factor= .07, min= 8):
    size = points * factor
    return max([size,min])

def fancy_dendrogram(linkage,names,name2color=False,threshold=False,self_thresh=False):

    # Make the dendrogram
    if threshold == False:
        scipy.cluster.hierarchy.dendrogram(linkage,labels=names,orientation='right')
    else:
        scipy.cluster.hierarchy.dendrogram(linkage,labels=names, color_threshold=threshold,\
                                            orientation='right')

    # Color the names
    if name2color != False:
        ax = plt.gca()
        xlbls = ax.get_ymajorticklabels()
        for lbl in xlbls:
            lbl.set_color(name2color[lbl.get_text()])

    # Add the threshold
    if threshold:
        plt.axvline(x=threshold, c='k')
    if self_thresh:
        plt.axvline(x=self_thresh, c='red', linestyle='dotted', lw=1)

    g = plt.gcf()
    return g

def normalize(df):
    result = df.copy()
    for feature_name in df.columns:
        if feature_name == 'genome':
            continue
        max_value = max(df[feature_name].max(),0)
        result[feature_name] = [max((x / max_value) if max_value != 0 else 0,0) for x in result[feature_name].tolist()]

    return result

def gen_color_list(names,name2cluster):
    '''
    Make a list of colors the same length as names, based on their cluster
    '''
    cm = plt.get_cmap('gist_rainbow')

    # 1. generate cluster to color
    cluster2color = {}
    clusters = set(name2cluster.values())
    NUM_COLORS = len(clusters)
    for cluster in clusters:
        try:
            cluster2color[cluster] = cm(1.*int(cluster)/NUM_COLORS)
        except:
            cluster2color[cluster] = cm(1.*int(str(cluster).split('_')[1])/NUM_COLORS)

    #2. generate list of colors
    colors = []
    for name in names:
        colors.append(cluster2color[name2cluster[name]])

    return colors

def gen_color_dictionary(names,name2cluster):
    '''
    Make the dictionary name2color
    '''
    cm = rand_cmap(len(set(name2cluster.values())),type='soft')

    # 1. generate cluster to color
    cluster2color = {}
    clusters = set(name2cluster.values())
    NUM_COLORS = len(clusters)
    for cluster in clusters:
        try:
            cluster2color[cluster] = cm(1.*int(cluster)/NUM_COLORS)
        except:
            cluster2color[cluster] = cm(1.*int(str(cluster).split('_')[1])/NUM_COLORS)

    #2. name to color
    name2color = {}
    for name in names:
        name2color[name] = cluster2color[name2cluster[name]]

    return name2color

def gen_colors(name2cluster):
    name2color = {}
    NUM_COLORS = len(name2cluster)
    cm = plt.get_cmap('gist_rainbow')

    for i,tax in enumerate(taxa):
        t2c[tax] = cm(1.*i/NUM_COLORS)

    return t2c

def comp_cluster(c):
    first = int(c.split('_')[0])
    dec = float(c.split('_')[1])
    return first + dec/100

def rand_cmap(nlabels, type='bright', first_color_black=True, last_color_black=False, verbose=False):
    """
    Creates a random colormap to be used together with matplotlib. Useful for segmentation tasks
    :param nlabels: Number of labels (size of colormap)
    :param type: 'bright' for strong colors, 'soft' for pastel colors
    :param first_color_black: Option to use first color as black, True or False
    :param last_color_black: Option to use last color as black, True or False
    :param verbose: Prints the number of labels and shows the colormap. True or False
    :return: colormap for matplotlib
    """
    from matplotlib.colors import LinearSegmentedColormap
    import colorsys
    import numpy as np


    if type not in ('bright', 'soft'):
        print ('Please choose "bright" or "soft" for type')
        return

    if verbose:
        print('Number of labels: ' + str(nlabels))

    # Generate color map for bright colors, based on hsv
    if type == 'bright':
        randHSVcolors = [(np.random.uniform(low=0.0, high=1),
                          np.random.uniform(low=0.2, high=1),
                          np.random.uniform(low=0.9, high=1)) for i in range(nlabels)]

        # Convert HSV list to RGB
        randRGBcolors = []
        for HSVcolor in randHSVcolors:
            randRGBcolors.append(colorsys.hsv_to_rgb(HSVcolor[0], HSVcolor[1], HSVcolor[2]))

        if first_color_black:
            randRGBcolors[0] = [0, 0, 0]

        if last_color_black:
            randRGBcolors[-1] = [0, 0, 0]

        random_colormap = LinearSegmentedColormap.from_list('new_map', randRGBcolors, N=nlabels)

    # Generate soft pastel colors, by limiting the RGB spectrum
    if type == 'soft':
        low = 0.6
        high = 0.95
        randRGBcolors = [(np.random.uniform(low=low, high=high),
                          np.random.uniform(low=low, high=high),
                          np.random.uniform(low=low, high=high)) for i in range(nlabels)]

        if first_color_black:
            randRGBcolors[0] = [0, 0, 0]

        if last_color_black:
            randRGBcolors[-1] = [0, 0, 0]
        random_colormap = LinearSegmentedColormap.from_list('new_map', randRGBcolors, N=nlabels)

    # Display colorbar
    if verbose:
        from matplotlib import colors, colorbar
        from matplotlib import pyplot as plt
        fig, ax = plt.subplots(1, 1, figsize=(15, 0.5))

        bounds = np.linspace(0, nlabels, nlabels + 1)
        norm = colors.BoundaryNorm(bounds, nlabels)

        cb = colorbar.ColorbarBase(ax, cmap=random_colormap, norm=norm, spacing='proportional', ticks=None,
                                   boundaries=bounds, format='%1i', orientation=u'horizontal')

    return random_colormap

def test_clustering():
    print("You should make some test cases here!")

if __name__ == '__main__':
	test_clustering()
