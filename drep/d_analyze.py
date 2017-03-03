#!/usr/bin/env python3

import matplotlib
matplotlib.use('Agg')

import pandas as pd
import os
import seaborn as sns
from matplotlib import pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
#matplotlib.rcParams['pdf.fonttype'] = 42
#matplotlib.rcParams['ps.fonttype'] = 42
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
        logging.debug("calling cluster_vis_wrapper")

        # If taxonomy info exists, add it
        if kwargs.get('include_taxonomy',True):
            Bdb = wd.get_db('Bdb')
            if 'taxonomy' in Bdb:
                genome2taxonomy = Bdb.set_index('genome')['taxonomy'].to_dict()
                kwargs['genome2taxonomy'] = genome2taxonomy

        cluster_vis_wrapper(wd, **kwargs)

    if kwargs.get('cluster') != None:
        logging.debug("calling cluster_test_wrapper")
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
    #print("making plots {0}".format(' '.join(sorted(to_plot))))

    # 1) Primary clustering dendrogram
    if '1' in to_plot:
        # Load the required data
        Mdb = wd.get_db('Mdb')
        Cdb = wd.get_db('Cdb')
        Pcluster = wd.get_primary_linkage()
        Plinkage = Pcluster['linkage']

        clust_args = wd.arguments['cluster']
        PL_thresh = clust_args.get('P_ani', False)
        if PL_thresh != False:
            PL_thresh = 1-PL_thresh

        # Make the plot
        logging.info("Plotting primary dendrogram...")
        plot_MASH_dendrogram(Mdb, Cdb, Plinkage, threshold = PL_thresh,\
                        plot_dir = plot_dir)

    # 2) Secondary clustering dendrogram
    if '2' in to_plot:
        logging.info("Plotting secondary dendrograms...")

        if 'Blank' in wd.get_db('Ndb'):
            logging.error("Nevermind- you don't have secondary clusters. Skipping plot 2")
        else:
            plot_secondary_dendrograms(wd, plot_dir, **kwargs)

    # 3) Secondary clusters heatmap
    if '3' in to_plot:
        # Load the required data
        logging.info("Plotting secondary clusters heatmaps...")
        logging.error("This part hasn't been written yet ¯\_(ツ)_/¯")

    # 4) Comparison scatterplots
    if '4' in to_plot:
        sns.set_style('whitegrid')

        # Load the required data
        Ndb = wd.get_db('Ndb')
        Mdb = wd.get_db('Mdb')
        Cdb = wd.get_db('Cdb')

        # Make the plot
        logging.info("Plotting Scatterplots...")

        if 'Blank' in wd.get_db('Ndb'):
            logging.error("Nope- you don't have secondary clusters. Skipping")
        else:
            plot_scatterplots(Mdb, Ndb, Cdb, plot_dir = plot_dir)

    # 5) Complex bin scorring
    if '5' in to_plot:
        # Load the required data
        Sdb = wd.get_db('Sdb')
        Cdb = wd.get_db('Cdb')
        Wdb = wd.get_db('Wdb')
        Chdb = wd.get_db('Chdb')

        # Make the plot
        logging.info("Plotting bin scorring plot...")
        plot_winner_scoring_complex(Wdb, Sdb, Cdb, Chdb, plot_dir = plot_dir, **kwargs)

    # 6) Winning plot
    if '6' in to_plot:
        # Load the required data
        Wdb = wd.get_db('Wdb')
        Chdb = wd.get_db('Chdb')
        Wndb = wd.get_db('Wndb')
        Wmdb = wd.get_db('Wmdb')
        Widb = wd.get_db('Widb')

        # Make the plot
        logging.info("Plotting winning genomes plot...")
        plot_winners(Wdb, Chdb, Wndb, Wmdb, Widb, plot_dir = plot_dir, **kwargs)

def cluster_test_wrapper(wd, **kwargs):

    # Validate arguments
    cluster = kwargs.get('cluster')
    comp_method = kwargs.get('clustering_method','ANIn')
    assert comp_method in ['ANIn','gANI']
    clust_method = kwargs.get('clusterAlg')
    threshold = kwargs.pop('threshold',None)
    cov_thresh = float(kwargs.get('minimum_coverage'))
    if threshold != None: threshold = 1- float(threshold)

    # Make a bdb listing the genomes to cluster
    Cdb = wd.get_db('Cdb')
    Bdb = wd.get_db('Bdb')
    genomes = Cdb['genome'][Cdb['primary_cluster'] == int(cluster)].tolist()
    bdb = Bdb[Bdb['genome'].isin(genomes)]

    # Get taxonomy if applicable
    if 'taxonomy' in Bdb:
        genome2taxonomy = Bdb.set_index('genome')['taxonomy'].to_dict()
        kwargs['genome2taxonomy'] = genome2taxonomy

    # Make the comparison database
    Xdb = dClust.compare_genomes(bdb,comp_method,wd,**kwargs)

    # Remove values without enough coverage
    if comp_method == 'ANIn':
        Xdb.loc[Xdb['alignment_coverage'] <= cov_thresh, 'ani'] = 0

    # Make it symmetrical
    Xdb['av_ani'] = Xdb.apply(lambda row: dClust.average_ani (row,Xdb),axis=1)
    Xdb['dist'] = 1 - Xdb['av_ani']
    db = Xdb.pivot("reference","querry","dist")

    # Cluster it
    if threshold == None:
        threshold = float(0)
    cdb, linkage = dClust.cluster_hierarchical(db, linkage_method = clust_method, \
                            linkage_cutoff = threshold)

    # Make the plot
    names = list(db.columns)
    if comp_method == 'ANIn':
        kwargs['self_thresh'] = get_highest_self(Xdb, names)
    kwargs['threshold'] = threshold
    kwargs['title_string']="Primary_cluster_{0}_{1}".format(cluster,clust_method)
    kwargs['subtitle_string'] = "Comp method: {0}    ".format(comp_method) +\
            "Clust method: {0}    Min cov: {1}".format(clust_method, cov_thresh)
    kwargs['name2cluster'] = cdb.set_index('genome')['cluster'].to_dict()


    plot_clustertest(linkage, names, wd, **kwargs)


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
                        logging.error("Can't interpret argument {0}- quitting".format(arg))
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

def _make_heatmap(db):
    g = sns.heatmap(db)
    labels = list(db.columns)

    # Adjust figure size
    fig = plt.gcf()
    fig.set_size_inches(x_fig_size(len(labels), factor=.5),x_fig_size(len(labels), factor=.5))
    plt.subplots_adjust(left=0.3)
    plt.subplots_adjust(bottom=0.3)

"""
SCATTER PLOTS
"""

def plot_scatterplots(Mdb, Ndb, Cdb, plot_dir=False):
    save = False

    # Initialize a .pdf
    if plot_dir != False:
        pp = PdfPages(plot_dir + 'Clustering_scatterplots.pdf')
        save = True

    g = plot_MASH_vs_ANIn_ani(Mdb,Ndb,exclude_zero_MASH=False)
    if save: pp.savefig(g)
    plt.show()

    g = plot_MASH_vs_secondary_ani(Mdb,Ndb,Cdb,exclude_zero_MASH=False)
    if save: pp.savefig(g)
    plt.show()

    g = plot_MASH_vs_ANIn_cov(Mdb,Ndb)
    if save: pp.savefig(g)
    plt.show()

    g = plot_ANIn_vs_ANIn_cov(Ndb)
    if save: pp.savefig(g)
    plt.show()

    try:
        g = plot_MASH_vs_len(Mdb,Ndb)
        if save: pp.savefig(g)
        plt.show()
    except:
        pass

    try:
        g = plot_ANIn_vs_len(Mdb,Ndb)
        if save: pp.savefig(g)
        plt.show()
    except:
        pass

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
    plt.gcf().suptitle('MASH vs ANIn comparisons (all)')
    plt.subplots_adjust(top=0.9)
    #plt.subplots_adjust(left=0.2)
    return plt.gcf()

def plot_MASH_vs_secondary_ani(Mdb,Ndb,Cdb,exclude_zero_MASH=True):
    plt.close('all')
    Xdb = pd.DataFrame()

    # Make a database of all comparisons
    mdb = Mdb.copy()
    mdb.rename(columns={'genome1':'querry','genome2':'reference',
                        'similarity':'MASH_ANI'},inplace=True)
    if exclude_zero_MASH:
        mdb= mdb[mdb['MASH_ANI'] > 0]
    db = pd.merge(mdb,Ndb)
    db.rename(columns={'ani':'ANIn'},inplace=True)

    # Filter to only comparisons within secondary clusters
    g2c = Cdb.set_index('genome')['secondary_cluster'].to_dict()
    db['ref_secondary_cluster'] = db['reference'].map(g2c)
    db['qu_secondary_cluster'] = db['querry'].map(g2c)
    db = db[db['ref_secondary_cluster'] == db['qu_secondary_cluster']]

    g = sns.jointplot(x='ANIn',y='MASH_ANI',data=db)
    plt.gcf().suptitle('MASH vs ANIn comparisons (within secondary clusters only)')
    plt.subplots_adjust(top=0.9)
    #plt.subplots_adjust(left=0.2)
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

    # Make a decending xaxis
    axs = plt.gcf().get_axes()
    plt.sca(axs[0])
    plt.xlim(db['MASH_ANI'].max(),db['MASH_ANI'].min())

    plt.gcf().suptitle('MASH vs length difference of genomes compared')
    plt.subplots_adjust(top=0.9)
    plt.subplots_adjust(left=0.2)

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

    # Make a decending xaxis
    axs = plt.gcf().get_axes()
    plt.sca(axs[0])
    plt.xlim(db['ANIn'].max(),db['ANIn'].min())

    plt.gcf().suptitle('ANIn vs length difference of genomes compared')
    plt.subplots_adjust(top=0.9)
    plt.subplots_adjust(left=0.2)
    return plt.gcf()

"""
CLUSETER PLOTS
"""

def plot_MASH_dendrogram(Mdb, Cdb, linkage, threshold = False, plot_dir = False, **kwargs):
    sns.set_style('white',{'axes.grid': False})

    db = Mdb.pivot("genome1","genome2","similarity")
    names = list(db.columns)
    name2cluster = Cdb.set_index('genome')['primary_cluster'].to_dict()
    name2color = gen_color_dictionary(names, name2cluster)

    # Make the dendrogram
    g = fancy_dendrogram(linkage,names,name2color,threshold=threshold)
    plt.title('MASH clustering')
    plt.xlabel('MASH Average Nucleotide Identity (ANI)')
    #plt.xlim([0,.4])

    sns.despine(left=True,top=True,right=True,bottom=False)

    # Adjust the figure size
    fig = plt.gcf()
    fig.set_size_inches(10,x_fig_size(len(names),factor=.2))
    plt.subplots_adjust(left=0.3)

    # Adjust the x labels
    plt.tick_params(axis='both', which='major', labelsize=8)
    axes = plt.gca()
    labels = axes.xaxis.get_majorticklocs()
    for i, label in enumerate(labels):
        labels[i] = (1 - float(label)) * 100
    axes.set_xticklabels(labels)

    # Add cluster to the y axis
    g2c = Cdb.set_index('genome')['secondary_cluster'].to_dict()
    axes = plt.gca()
    labels = [item.get_text() for item in axes.get_yticklabels()]
    for i, label in enumerate(labels):
        labels[i] = "{0} ({1})".format(label, g2c[label])
    axes.set_yticklabels(labels)

    # Save the figure
    if plot_dir != None:
        plt.savefig(plot_dir + 'Primary_clustering_dendrogram.pdf',format="pdf",\
            transparent=True, bbox_inches='tight')
    plt.show()
    plt.close('all')

def plot_secondary_dendrograms(wd, plot_dir, **kwargs):
    save = False
    win = False

    # Initialize a .pdf
    if plot_dir != False:
        pp = PdfPages(plot_dir + 'Secondary_clustering_dendrograms.pdf')
        save = True
    else:
        save = False

    # Load winner database if it exists
    if wd.hasDb('Wdb'):
        Wdb = wd.get_db('Wdb')
        winners = Wdb['genome'].unique()
        kwargs['winners'] = winners

    # Load secondary databases if they exist (Gdb, Ndb)
    Sdb = {}
    if wd.hasDb('Ndb'):
        Sdb['ANIn'] = wd.get_db('Ndb')
    if wd.hasDb('Gdb'):
        Sdb['gANI'] = wd.get_db('Gdb')

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
        min_cov = args['minimum_coverage']

        kwargs['threshold'] = threshold
        kwargs['title_string'] = 'Primary cluster {0}'.format(cluster)
        kwargs['subtitle_string'] = "Comp method: {0}    ".format(alg) +\
                "Clust method: {0}    Min cov: {1}".format(clust_alg, min_cov)

        # Get name2cluster
        names = list(db.columns)
        name2cluster = Cdb.set_index('genome')['secondary_cluster'].to_dict()
        for name in names: # Handle the case where you deleted a secondary cluster
            if name not in Cdb['genome'].tolist():
                name2cluster[name] = '0_0'
        kwargs['name2cluster'] = name2cluster

        # Get the highest self-comparison
        #if alg in Sdb:
            #kwargs['self_thresh'] = get_highest_self(Sdb[alg], names)
        kwargs['self_thresh'] = get_highest_self(wd.get_db('Ndb'), names)

        # Make the dendrogram
        _make_special_dendrogram(linkage,names,**kwargs)

        # Save the file
        fig = plt.gcf()
        if save == True:
            pp.savefig(fig)
        plt.show()
        plt.close(fig)

    pp.close()
    plt.close('all')

'''
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
            plt.xlim([0,3*threshold])

        # Adjust the figure size
        fig = plt.gcf()
        fig.set_size_inches(10,x_fig_size(len(names),factor=.5))
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
'''


def plot_clustertest(linkage, names, wd, **kwargs):
    '''
    names can be gotten like:
    db = db.pivot("reference","querry","ani")
    names = list(db.columns)
    '''

    # Make the plot directory
    plot_dir = wd.location + '/figures/cluster_tests/'
    if not os.path.exists(plot_dir):
        os.makedirs(plot_dir)

    # Make the dendrogram
    _make_special_dendrogram(linkage,names,**kwargs)

    # Save the dendrogram
    fig = plt.gcf()
    plt.savefig("{0}{1}.pdf".format(plot_dir, kwargs['title_string']))
    plt.show()
    plt.close(fig)

'''
names can be gotten like:
db = db.pivot("reference","querry","ani")
names = list(db.columns)
'''
def _make_dendrogram(linkage, names, **kwargs):
    threshold = kwargs.get('threshold',False)
    title = kwargs.get('title',False)

    # Make the dendrogram
    g = fancy_dendrogram(linkage,names,threshold=threshold)
    plt.title(title)
    plt.xlabel(kwargs.get('xlabel','ANI'))
    if threshold != False:
        plt.xlim([0,2*threshold])

    # Adjust the figure size
    fig = plt.gcf()
    fig.set_size_inches(10,x_fig_size(len(names),factor=.2))
    plt.subplots_adjust(left=0.3)

    # Adjust the x labels
    plt.tick_params(axis='both', which='major', labelsize=8)
    axes = plt.gca()
    labels = axes.xaxis.get_majorticklocs()
    for i, label in enumerate(labels):
        labels[i] = (1 - float(label)) * 100
    axes.set_xticklabels(labels)

'''
names can be gotten like:
db = db.pivot("reference","querry","ani")
names = list(db.columns)
'''
def _make_special_dendrogram(linkage, names, **kwargs):
    # Load possible kwargs
    name2cluster = kwargs.get('name2cluster',False)
    self_thresh = kwargs.get('self_thresh',False)
    threshold = kwargs.get('threshold',False)
    title_string = kwargs.get('title_string','')
    subtitle_string = kwargs.get('subtitle_string','')
    winners = kwargs.get('winners',False)
    genome2taxonomy = kwargs.get('genome2taxonomy',False)

    # Make special things
    if name2cluster != False:
        name2color = gen_color_dictionary(names, name2cluster)
    else:
        name2color = False

    # Make the dendrogram
    sns.set_style('whitegrid')
    g = fancy_dendrogram(linkage,names,name2color,threshold=threshold,self_thresh =\
                        self_thresh)

    # Add the title and subtitle
    plt.suptitle(title_string, y=1, fontsize=18)
    plt.title(subtitle_string, fontsize=10)

    # Adjust the x-axis
    plt.xlabel('Average Nucleotide Identity (ANI)')
    if threshold != False:
        plt.xlim([0,3*threshold])
    plt.tick_params(axis='both', which='major', labelsize=12)
    axes = plt.gca()
    labels = axes.xaxis.get_majorticklocs()
    for i, label in enumerate(labels):
        labels[i] = (1 - float(label)) * 100
    axes.set_xticklabels(labels)
    plt.gca().yaxis.grid(False)

    # Adjust the figure size
    fig = plt.gcf()
    fig.set_size_inches(10,x_fig_size(len(names),factor=.5))
    plt.subplots_adjust(left=0.5)

    # Mark winning ones
    if type(winners) is not bool:
        ax = plt.gca()
        labels = [item.get_text() for item in ax.get_yticklabels()]
        for i, label in enumerate(labels):
            if label in winners: labels[i] = label + ' *'
        ax.set_yticklabels(labels)

    # Add taxonomy
    if genome2taxonomy != False:
        g2t = genome2taxonomy
        axes = plt.gca()
        labels = [item.get_text() for item in axes.get_yticklabels()]
        for i, label in enumerate(labels):
            labels[i] = "{0}\n{1}".format(label, g2t[label.replace(' *','')])
        axes.set_yticklabels(labels)

"""
WINNER PLOTS
"""

def plot_winner_scoring_complex(Wdb, Sdb, Cdb, Chdb, plot_dir= False, **kwargs):
    save = False

    # Initialize a .pdf
    if plot_dir != False:
        pp = PdfPages(plot_dir + 'Cluster_scoring.pdf')
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

'''
db is the database to plot- must contain 'genome' and all columns listed in 'bars'
bars is all of the columns in the database to become bars
for taxonomy, put genome2taxonomy in kwargs
'''
def _make_scoring_plot(db, bars,**kwargs):
    sns.set_style('whitegrid')

    # Make the normalized bar plot
    nd = normalize(db)
    d = pd.melt(nd, id_vars=['genome'], value_vars=bars)
    g = sns.barplot(data=d, y='genome', x='value', hue='variable')

    # Get a list of the un-normalized values
    x = pd.melt(db, id_vars=['genome'], value_vars=bars)
    vals = []
    for variable in x['variable'].unique():
        vals += [v for v in x['value'][x['variable'] == variable].tolist()]

    # Add un-normalized values to barplots
    i = 0
    for p in g.patches:
        g.annotate("{0:.1f}".format(vals[i]), (p.get_width(), p.get_y()+(p.get_height()/1.1) ), fontsize=8)
        i += 1

    # Add taxonomy if available
    axes = plt.gca()
    labels = [item.get_text() for item in axes.get_yticklabels()]
    if kwargs.get('genome2taxonomy',False) != False:
        g2t = kwargs.get('genome2taxonomy')
        for i, label in enumerate(labels):
            labels[i] = "{0}\n{1}".format(label, g2t[label.replace(' *','')])
        axes.set_yticklabels(labels)

    # Adjust labels
    plt.xlabel('Normalized Score')
    plt.legend(loc='lower right')
    plt.tick_params(axis='both', which='major', labelsize=8)

    # Adjust figure size
    fig = plt.gcf()
    fig.set_size_inches(12,x_fig_size(len(labels), factor=1))
    plt.subplots_adjust(left=0.5)

def plot_winners(Wdb, Chdb, Wndb, Wmdb, Widb, plot_dir= False, **kwargs):
    save = False

    # Initialize a .pdf
    if plot_dir != False:
        pp = PdfPages(plot_dir + 'Winning_genomes.pdf')
        save = True

    # Make piecharts
    labels = []
    sizes = []
    for com in Widb['completeness_metric'].unique():
        d = Widb[Widb['completeness_metric'] == com]
        labels.append(com)
        sizes.append(len(d['genome'].unique()))
    labels = annotate_labels(labels,'comp')
    _make_piechart(labels,sizes)
    plt.title('Overall Winner Completeness')

    # Save this page
    if save == True:
        fig = plt.gcf()
        pp.savefig(fig)
    plt.show()
    plt.close(fig)

    labels = []
    sizes = []
    for com in Widb['contamination_metric'].unique():
        d = Widb[Widb['contamination_metric'] == com]
        labels.append(com)
        sizes.append(len(d['genome'].unique()))
    labels = annotate_labels(labels,'con')
    _make_piechart(labels,sizes)
    plt.title('Overall Winner Contamination')

    # Save this page
    if save == True:
        fig = plt.gcf()
        pp.savefig(fig)
    plt.show()
    plt.close(fig)

    # Figure out what you're going to show
    bars = kwargs.get('to_show',['score','N50 (scaffolds)','Completeness',\
            'Contamination','Genome size (bp)', 'Strain heterogeneity'])

    # Make a db for the winners
    Chdb['genome'] = Chdb['Bin Id']
    d = Wdb.sort_values('score', ascending=False)
    d = d.merge(Chdb, how='left', on= 'genome')
    d = d[bars + ['genome']]

    # Make the scoring plot
    _make_scoring_plot(d,bars,**kwargs)
    plt.title('Scoring of winning genomes')

    # Save this page
    if save == True:
        fig = plt.gcf()
        pp.savefig(fig)
    plt.show()
    plt.close(fig)

    if Wmdb is not None:
        # Make a MASH linkage for the dendrogram
        db = Wmdb.copy()
        db['dist'] = 1 - db['similarity']
        linkage_db = db.pivot("genome1","genome2","dist")
        names = list(linkage_db.columns)
        Cdb, linkage = dClust.cluster_hierarchical(linkage_db, linkage_method= 'average', \
                                    linkage_cutoff= 0)

        # Make the MASH dendrogram
        _make_dendrogram(linkage,names)
        plt.title('MASH dendrogram')

        # Save this page
        if save == True:
            fig = plt.gcf()
            pp.savefig(fig)
        plt.show()
        plt.close(fig)

    if Wndb is not None:
        # Make a ANIn linkage for the dendrogram
        d = Wndb.copy()
        d['av_ani'] = d.apply(lambda row: dClust.average_ani (row,d),axis=1)
        d['dist'] = 1 - d['av_ani']
        db = d.pivot("reference", "querry", "dist")
        names = list(db.columns)
        Cdb, linkage = dClust.cluster_hierarchical(db, linkage_method= 'average', \
                                    linkage_cutoff= 0)

        # Make the ANIn dendrogram
        _make_dendrogram(linkage,names)
        plt.title('ANIn dendrogram (NOT filtered for alignment length)')

        # Save this page
        if save == True:
            fig = plt.gcf()
            pp.savefig(fig)
        plt.show()
        plt.close(fig)

        # Make a ANIn linkage for the filtered dendrogram
        d = Wndb.copy()
        d.loc[d['alignment_coverage'] <= 0.1, 'ani'] = 0
        d['av_ani'] = d.apply(lambda row: dClust.average_ani (row,d),axis=1)
        d['dist'] = 1 - d['av_ani']
        db = d.pivot("reference", "querry", "dist")
        names = list(db.columns)
        Cdb, linkage = dClust.cluster_hierarchical(db, linkage_method= 'average', \
                                    linkage_cutoff= 0)

        # Make the ANIn dendrogram
        _make_dendrogram(linkage,names)
        plt.title('ANIn dendrogram (filtered for 10% alignment)')

        # Save this page
        if save == True:
            fig = plt.gcf()
            pp.savefig(fig)
        plt.show()
        plt.close(fig)

    # Save the .pdf
    if save:
        pp.close()
    plt.close('all')

"""
OTHER
"""

def get_highest_self(db, genomes, min = 1.0e-4):
    d = db[db['reference'].isin(genomes)]
    self_thresh = 1 - d['ani'][d['reference'] == d['querry']].min()

    # Because 0s don't show up on the graph
    if self_thresh == float(0):
        self_thresh =  min
    return self_thresh

def _make_piechart(labels,sizes):
    plt.pie(sizes,labels=labels,startangle=45,\
            autopct=make_autopct(sizes),\
            shadow = True)
    plt.axis('equal')

def make_autopct(values):
    def my_autopct(pct):
        total = sum(values)
        val = int(round(pct*total/100.0))
        return '{p:.2f}%  ({v:d})'.format(p=pct,v=val)
    return my_autopct

def annotate_labels(labels,how):
    if how == 'comp':
        labs = []
        for label in labels:
            if label == 'near':
                labs.append('near (>90%)')
            if label == 'perfect':
                labs.append('perfect (100%)')
            if label == 'substantial':
                labs.append('substantial (>70%)')
            if label == 'moderate':
                labs.append('moderate (>50%)')
            if label == 'partial':
                labs.append('partial (<50%)')
        return labs

    if how == 'con':
        labs = []
        for label in labels:
            if label == 'low':
                labs.append('low (<5%)')
            if label == 'none':
                labs.append('none (0%)')
            if label == 'medium':
                labs.append('medium (<10%)')
            if label == 'high':
                labs.append('high (<15%)')
            if label == 'very high':
                labs.append('very high (>15%)')
        return labs

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
        plt.axvline(x=threshold, c='k', linestyle='dotted')
    if self_thresh != False:
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
    cm = rand_cmap(len(set(name2cluster.values()))+1,type='bright')

    # 1. generate cluster to color
    cluster2color = {}
    clusters = set(name2cluster.values())
    NUM_COLORS = len(clusters)
    for cluster in clusters:
        try:
            #x = cm(1.*int(cluster)/NUM_COLORS)
            #print(x)
            cluster2color[cluster] = cm(1.*int(cluster)/NUM_COLORS)
        except:
            #x = cm(1.*int(cluster)/NUM_COLORS)
            #print(x)
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
