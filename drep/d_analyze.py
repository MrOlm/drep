#!/usr/bin/env python3
'''
d_analyze - a subset of drep

Make plots based on de-replication
'''
import matplotlib
matplotlib.use('Agg')

import logging
import math
import os

import pandas as pd
import seaborn as sns
import scipy.cluster.hierarchy
from sklearn import manifold
import numpy as np

from matplotlib import pyplot as plt
import matplotlib.ticker as ticker
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.patches as mpatches

import drep
import drep.d_cluster
import drep.d_filter

import traceback

#warnings.filterwarnings("ignore", category=MatplotlibDeprecationWarning)

def d_analyze_wrapper(wd, **kwargs):
    '''
    Controller for the dRep analyze operation

    Args:
        wd: The current workDirectory
        **kwargs: Command line arguments

    Keyword Args:
        plots: List of plots to make [list of ints, 1-6]

    Returns:
        Makes some plots
    '''

    # Load the workDirectory
    wd = drep.WorkDirectory.WorkDirectory(wd)

    # Figure out what plots to make
    options = ['1','2','3','4','5','6']
    to_plot = _parse_plot_options(options, kwargs.get('plots', None))
    logging.info("making plots {0}".format(', '.join(to_plot)))

    # Get the plot directory
    plot_dir = wd.get_dir('figures')

    # 1) Primary clustering dendrogram
    if '1' in to_plot:
        try:
            mash_dendrogram_from_wd(wd, plot_dir=plot_dir)
        except BaseException as e:
            logging.info('Failed to make plot #1: ' + str(e))
            traceback.print_exc()

    # 2) Secondary clustering dendrogram
    if '2' in to_plot:
        try:
            plot_secondary_dendrograms_from_wd(wd, plot_dir, **kwargs)
        except BaseException as e:
            logging.info('Failed to make plot #2: ' + str(e))
            traceback.print_exc()

    # 3) Secondary clusters MDS
    if '3' in to_plot:
        try:
            plot_secondary_mds_from_wd(wd, plot_dir, **kwargs)
        except BaseException as e:
            logging.info('Failed to make plot #3: ' + str(e))
            traceback.print_exc()

    # 4) Comparison scatterplots
    if '4' in to_plot:
        try:
            plot_scatterplots_from_wd(wd, plot_dir, **kwargs)
        except BaseException as e:
            logging.info('Failed to make plot #4: ' + str(e))
            traceback.print_exc()

    # 5) Complex bin scorring
    if '5' in to_plot:
        try:
            plot_binscoring_from_wd(wd, plot_dir, **kwargs)
        except BaseException as e:
            logging.info('Failed to make plot #5: ' + str(e))
            traceback.print_exc()

    # 6) Winning plot
    if '6' in to_plot:
        try:
            plot_winners_from_wd(wd, plot_dir, **kwargs)
        except BaseException as e:
            logging.info('Failed to make plot #6: ' + str(e))
            traceback.print_exc()


def mash_dendrogram_from_wd(wd, plot_dir=False):
    '''
    From the wd and kwargs, call plot_MASH_dendrogram

    Args:
        wd: WorkDirectory
        plot_dir (optional): Location to store figure

    Returns:
        Shows plot, makes a plot in the plot_dir
    '''
    # Load the required data
    try:
        Mdb = wd.get_db('Mdb', return_none=False, forPlotting=True)
        Cdb = wd.get_db('Cdb', return_none=False)
        Pcluster = wd.get_primary_linkage()
        Plinkage = Pcluster['linkage']
        clust_args = wd.arguments['cluster']
        PL_thresh = clust_args.get('P_ani', False)
        if PL_thresh != False:
            PL_thresh = 1-PL_thresh
    except:
        logging.error("Skipping plot 1 - you don't have all required dataframes")
        return

    # Make the plot
    logging.info("Plotting primary dendrogram")
    plot_MASH_dendrogram(Mdb, Cdb, Plinkage, threshold = PL_thresh,\
                    plot_dir = plot_dir)

def plot_secondary_dendrograms_from_wd(wd, plot_dir, **kwargs):
    '''
    From the wd and kwargs, make the secondary dendrograms

    Args:
        wd: WorkDirectory
        plot_dir (optional): Location to store figure

    Returns:
        Makes plot
    '''

    # Initialize a .pdf
    if plot_dir != False:
        pp = PdfPages(plot_dir + 'Secondary_clustering_dendrograms.pdf')
        save = True
    else:
        save = False

    # Load required databases
    try:
        Ndb = wd.get_db('Ndb', return_none=False)
        Cdb = wd.get_db('Cdb', return_none=False)
    except:
        logging.error("Skipping plot 2 - you don't have all required dataframes")
        return
    logging.info("Plotting secondary dendrograms")

    # Load winner database if it exists
    if wd.hasDb('Wdb'):
        Wdb = wd.get_db('Wdb')
        winners = Wdb['genome'].unique()
        kwargs['winners'] = winners

    # Load genome 2 taxonomy if it exists
    genome2taxonomy = _get_genome2taxonomy(wd)
    kwargs['genome2taxonomy'] = genome2taxonomy

    # For every cluster:
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

        kwargs['self_thresh'] = get_highest_self(Ndb, names)

        # Make the dendrogram
        _make_special_dendrogram(linkage, names, **kwargs)

        # Save the file
        fig = plt.gcf()
        if save == True:
            pp.savefig(fig)
        plt.show()
        plt.close(fig)

    pp.close()
    plt.close('all')

def plot_secondary_mds_from_wd(wd, plot_dir, **kwargs):
    '''
    Make a .pdf of MDS of each cluster

    Args:
        wd: WorkDirectory
        plot_dir (optional): Location to store figure

    Returns:
        Makes plot
    '''
    # Load required databases
    try:
        Ndb = wd.get_db('Ndb', return_none=False)
        Cdb = wd.get_db('Cdb', return_none=False)
    except:
        logging.error("Skipping plot 3 - you don't have all required dataframes")
        return

    logging.info("Plotting MDS plot")
    # initialize a .pdf
    if plot_dir != False:
        pp = PdfPages(plot_dir + 'Secondary_clustering_MDS.pdf')
        save = True
    else:
        save = False

    # for every cluster:
    Cdb = wd.get_db('Cdb')
    for cluster in sorted(Cdb['primary_cluster'].unique()):
        d = Cdb[Cdb['primary_cluster'] == cluster]

        # Skip if it's a singleton
        if len(d['genome'].unique()) == 1:
            continue

        # Load the linkage information
        linkI = wd.get_cluster("secondary_linkage_cluster_{0}".format(cluster))
        db = linkI['db']
        args = linkI['arguments']

        # Load name to cluster
        names = list(db.columns)
        name2cluster = Cdb.set_index('genome')['secondary_cluster'].to_dict()
        for name in names: # Handle the case where you deleted a secondary cluster
            if name not in Cdb['genome'].tolist():
                name2cluster[name] = '0_0'

        # Get the colors
        name2color = gen_color_dictionary(names, name2cluster)
        colors = [name2color[n] for n in names]

        # Make cluster 2 color
        cluster2color = {name2cluster[n]: name2color[n] for n in names}

        # make the mds plot
        _make_mds_plot(cluster, db, names, colors=colors, annotate=False,
                      cluster2color = cluster2color)

        # save the plot
        fig = plt.gcf()
        if save == True:
            pp.savefig(fig, bbox_inches='tight')
        plt.show()
        plt.close(fig)

    pp.close()
    plt.close('all')

def plot_scatterplots_from_wd(wd, plot_dir, **kwargs):
    '''
    From the wd and kwargs, call plot_scatterplots

    Args:
        wd: WorkDirectory
        plot_dir (optional): Location to store figure

    Returns:
        Shows plot, makes a plot in the plot_dir
    '''
    # Load the required data
    try:
        Ndb = wd.get_db('Ndb', return_none=False)
        Mdb = wd.get_db('Mdb', return_none=False)
        Cdb = wd.get_db('Cdb', return_none=False)
    except:
        logging.error("Skipping plot 4 - you don't have all required dataframes")
        return

    # Make the plot
    logging.info("Plotting scatterplots")
    plot_scatterplots(Mdb, Ndb, Cdb, plot_dir = plot_dir)

def plot_binscoring_from_wd(wd, plot_dir, **kwargs):
    '''
    From the wd and kwargs, call plot_winner_scoring_complex

    Args:
        wd: WorkDirectory
        plot_dir (optional): Location to store figure

    Returns:
        Shows plot, makes a plot in the plot_dir
    '''
    # Load the required data
    try:
        Sdb = wd.get_db('Sdb', return_none=False)
        Cdb = wd.get_db('Cdb', return_none=False)
        Wdb = wd.get_db('Wdb', return_none=False)
        Bdb = wd.get_db('Bdb', return_none=False)
    except:
        logging.error("Skipping plot 5 - you don't have all required dataframes")
        return

    # Deal with genome quality
    Gdb = drep.d_filter._get_run_genomeInfo(wd, Bdb, no_run=True)

    # Only keep things you want
    Gdb = Gdb[[c for c in Gdb.columns if c in ['genome', 'location', 'N50', 'length', 'completeness', 'contamination', 'strain_heterogeneity']]]

    # Make the plot
    logging.info("Plotting bin scorring plot")
    plot_winner_scoring_complex(Wdb, Sdb, Cdb, Gdb, plot_dir = plot_dir, **kwargs)

def plot_winners_from_wd(wd, plot_dir, **kwargs):
    '''
    From the wd and kwargs, call plot_winners

    Args:
        wd: WorkDirectory
        plot_dir: Location to store figure

    Returns:
        Shows plot, makes a plot in the plot_dir
    '''
    # Load the required data
    try:
        Wdb = wd.get_db('Wdb', return_none=False)
        Bdb = wd.get_db('Bdb', return_none=False)
    except:
        logging.error("Skipping plot 6 - you don't have all required dataframes")
        return

    # Deal with genome quality
    Gdb = drep.d_filter._get_run_genomeInfo(wd, Bdb, no_run=True)

    # Only keep things you want
    Gdb = Gdb[[c for c in Gdb.columns if c in ['genome', 'location', 'N50', 'length', 'completeness', 'contamination', 'strain_heterogeneity']]]

    # Get optional data
    Wndb = wd.get_db('Wndb')
    Wmdb = wd.get_db('Wmdb')
    Widb = wd.get_db('Widb')

    # Make the plot
    logging.info("Plotting winning genomes plot...")
    plot_winners(Wdb, Gdb, Wndb, Wmdb, Widb, plot_dir = plot_dir, **kwargs)

def plot_scatterplots(Mdb, Ndb, Cdb, plot_dir=False):
    '''
    Make scatterplots comparing genome comparison algorithms

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

    Args:
        Mdb: DataFrame of Mash comparison results
        Ndb: DataFrame of secondary clustering results
        Cdb: DataFrame of Clustering results
        plot_dir (optional): Location to store plot

    Return:
        Makes and shows plot
    '''
    sns.set_style('whitegrid')

    # Initialize a .pdf
    if plot_dir != False:
        pp = PdfPages(plot_dir + 'Clustering_scatterplots.pdf')
        save = True
    else:
        save = False

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

def plot_MASH_vs_ANIn_ani(Mdb, Ndb, exclude_zero_MASH=True):
    '''
    Makes plot and retuns plt.cgf()

    All parameters are obvious
    '''
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
    '''
    Makes plot and retuns plt.cgf()

    All parameters are obvious
    '''
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
    '''
    Makes plot and retuns plt.cgf()

    All parameters are obvious
    '''
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
    '''
    Makes plot and retuns plt.cgf()

    All parameters are obvious
    '''
    plt.close('all')
    db = Ndb.copy()
    db.rename(columns={'alignment_coverage':'ANIn_alignment_coverage','ani':'ANIn'},inplace=True)
    g = sns.jointplot(y='ANIn',x='ANIn_alignment_coverage',data=db)
    return plt.gcf()

def plot_MASH_vs_len(Mdb,Ndb,exclude_zero_MASH=True):
    '''
    Makes plot and retuns plt.cgf()

    All parameters are obvious
    '''
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
    '''
    Makes plot and retuns plt.cgf()

    All parameters are obvious
    '''
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

def plot_MASH_dendrogram(Mdb, Cdb, linkage, threshold=False, plot_dir=False):
    '''
    Make a dendrogram of the primary clustering

    Args:
        Mdb: DataFrame of Mash comparison results; make sure loaded not as categories
        Cdb: DataFrame of Clustering results
        linkage: Result of scipy.cluster.hierarchy.linkage
        threshold (optional): Line to plot on x-axis
        plot_dir (optional): Location to store plot

    Returns:
        Makes and shows plot
    '''
    sns.set_style('white',{'axes.grid': False})

    if Mdb['genome1'].dtype.name == 'category':
        logging.error("WARNING: Primary dendrogram labels may be shuffled! Load as csv to prevent this")

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
    fig.set_size_inches(10,_x_fig_size(len(names),factor=.2))
    plt.subplots_adjust(left=0.3)

    # Adjust the x labels
    plt.tick_params(axis='both', which='major', labelsize=8)
    axes = plt.gca()
    labels = axes.xaxis.get_majorticklocs()
    for i, label in enumerate(labels):
        labels[i] = float("{0:.2f}".format((1 - float(label)) * 100))
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
        plt.savefig(os.path.join(plot_dir, 'Primary_clustering_dendrogram.pdf'),\
            format="pdf", transparent=True, bbox_inches='tight')
    plt.show()
    plt.close('all')

"""
WINNER PLOTS
"""

def plot_winner_scoring_complex(Wdb, Sdb, Cdb, Gdb, plot_dir= False, **kwargs):
    '''
    Make a plot showing the genome scoring for all genomes

    Args:
        Wdb: DataFrame of winning dereplicated genomes
        Sdb: Scores of all genomes
        Cdb: DataFrame of Clustering results
        Gdb: DataFrame of genome scoring information
        plot_dir (optional): Location to store plot

    Returns:
        makes plot
    '''
    # Set style
    sns.reset_orig()

    # Initialize a .pdf
    if plot_dir != False:
        pp = PdfPages(os.path.join(plot_dir, 'Cluster_scoring.pdf'))
        save = True
    else:
        save = False

    # Figure out what you're going to show
    bars = _get_toshow(Gdb)
    bars += ['score']

    # Get winners
    winners = list(Wdb['genome'].unique())

    for cluster in sorted(Cdb['secondary_cluster'].unique(), key=lambda x: _comp_cluster(x)):
        # Make a db for this cluster
        d = Cdb[Cdb['secondary_cluster'] == cluster]
        d = d.merge(Sdb, how='left', on= 'genome')
        d = d.merge(Gdb, how='left', on= 'genome')
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
        plt.legend(loc=(0,0), fancybox=True, framealpha=0.5)
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
        fig.set_size_inches(12,_x_fig_size(len(labels), factor=1))
        plt.subplots_adjust(left=0.5)

        if save == True:
            pp.savefig(fig)
        plt.show()
        plt.close(fig)

    if save:
        pp.close()
    plt.close('all')

def plot_winners(Wdb, Gdb, Wndb, Wmdb, Widb, plot_dir= False, **kwargs):
    '''
    Make a bunch of plots about the de-replicated genomes

    THIS REALLY NEEDS IMPROVED UPON
    '''

    # Set style
    sns.reset_orig()

    # Initialize a .pdf
    if plot_dir != False:
        pp = PdfPages(plot_dir + 'Winning_genomes.pdf')
        save = True
    else:
        save = False

    # Make piecharts
    if Widb is not None:
        labels = []
        sizes = []
        for com in Widb['completeness_metric'].unique():
            d = Widb[Widb['completeness_metric'] == com]
            labels.append(com)
            sizes.append(len(d['genome'].unique()))
        labels = _annotate_labels(labels,'comp')

        if (len(labels) != len(sizes)) | (len(labels) == 0): # not sure when this would happen, but it does...
            logging.debug("len(labels) != len(sizes); {0} vs {1}".format(\
                len(labels), len(sizes)))

        else:
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
            labels = _annotate_labels(labels,'con')
            _make_piechart(labels,sizes)
            plt.title('Overall Winner Contamination')

            # Save this page
            if save == True:
                fig = plt.gcf()
                pp.savefig(fig)
            plt.show()
            plt.close(fig)

    # Figure out what you're going to show
    bars = _get_toshow(Gdb)
    bars += ['score']

    # Make a db for the winners
    d = Wdb.sort_values('score', ascending=False)
    d = d.merge(Gdb, how='left', on= 'genome')
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
        Cdb, linkage = drep.d_cluster.cluster_hierarchical(linkage_db, linkage_method= 'average', \
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
        drep.d_cluster.add_avani(d)
        #d['av_ani'] = d.apply(lambda row: drep.d_cluster.average_ani (row,d),axis=1)
        d['dist'] = 1 - d['av_ani']
        db = d.pivot("reference", "querry", "dist")
        names = list(db.columns)
        Cdb, linkage = drep.d_cluster.cluster_hierarchical(db, linkage_method= 'average', \
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
        drep.d_cluster.add_avani(d)
        #d['av_ani'] = d.apply(lambda row: drep.d_cluster.average_ani (row,d),axis=1)
        d['dist'] = 1 - d['av_ani']
        db = d.pivot("reference", "querry", "dist")
        names = list(db.columns)
        Cdb, linkage = drep.d_cluster.cluster_hierarchical(db, linkage_method= 'average', \
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

def calc_dist(x1, y1, x2, y2):
    '''
    Return distance from two points

    Args: self explainatory

    Returns:
        int: distance
    '''
    dist = math.hypot(x2 - x1, y2 - y1)
    return dist

def get_highest_self(db, genomes, min = 1.0e-4):
    '''
    Return the highest ANI value resulting from comparing a genome to itself
    '''
    d = db[db['reference'].isin(genomes)]
    self_thresh = 1 - d['ani'][d['reference'] == d['querry']].min()

    # Because 0s don't show up on the graph
    if self_thresh == float(0):
        self_thresh =  min
    return self_thresh

def _make_piechart(labels,sizes):
    '''
    Used by winner plot
    '''
    plt.pie(sizes,labels=labels,startangle=45,\
            autopct=_make_autopct(sizes),\
            shadow = True)
    plt.axis('equal')

def _make_autopct(values):
    '''
    Used by winner plot
    '''
    def my_autopct(pct):
        total = sum(values)
        val = int(round(pct*total/100.0))
        return '{p:.2f}%  ({v:d})'.format(p=pct,v=val)
    return my_autopct

def _annotate_labels(labels,how):
    '''
    Used by winner plot
    '''
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

def _x_fig_size(points, factor= .07, min= 8):
    '''
    Calculate how big the x of the figure should be
    '''
    size = points * factor
    return max([size,min])

def fancy_dendrogram(linkage,names,name2color=False,threshold=False,self_thresh=False):
    '''
    Make a fancy dendrogram
    '''
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
    '''
    Normalize all columns in df to 0-1 except 'genome' or 'location'

    Args:
        df: DataFrame

    Return:
        DataFrame: Nomralized
    '''
    result = df.copy()
    for feature_name in df.columns:
        if feature_name in ['genome', 'location']:
            continue
        if not np.issubdtype(df[feature_name].dtype, np.number):
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

def gen_color_dictionary(names, name2cluster):
    '''
    Make the dictionary name2color

    Args:
        names: key in the returned dictionary
        name2cluster: a dictionary of name to it's cluster

    Returns:
        dict: name -> color
    '''
    #cm = _rand_cmap(len(set(name2cluster.values()))+1,type='bright')
    vals = np.linspace(0,1,len(set(name2cluster.values()))+1)
    np.random.shuffle(vals)
    cm = plt.cm.colors.ListedColormap(plt.cm.jet(vals))

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

def _comp_cluster(c):
    '''
    Used in secondary dendrogram creation
    '''
    first = int(c.split('_')[0])
    dec = float(c.split('_')[1])
    return first + dec/100

def _rand_cmap(nlabels, type='bright', first_color_black=True, last_color_black=False, verbose=False):
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


def _parse_plot_options(options, args):
    '''
    Read user input and figure out a list of plots to make

    Args:
        options: list of possible plots to make (default [1-6])
        args: the command line passed in

    Returns:
        list: list of ints in the args
    '''
    to_plot = []

    if args[0] in ['all','a']:
        to_plot += options

    elif args == None:
        logging.error("No plots given!")
        sys.exit()
        return None

    else:
        for arg in args:
            if arg in options:
                to_plot.append(arg)
            else:
                for letter in arg:
                    if letter in options:
                        to_plot.append(letter)
                    else:
                        logging.error("Can't interpret plotting argument {0}! quitting".format(arg))
                        sys.exit()
    return to_plot

def _get_genome2taxonomy(wd):
    '''
    Return dictionary: genome -> taxonomy

    All based on Bdb at the moment

    Return False if can't do it
    '''
    try:
        Bdb = wd.get_db('Bdb')
        if 'taxonomy' in Bdb:
            genome2taxonomy = Bdb.set_index('genome')['taxonomy'].to_dict()
        return genome2taxonomy
    except:
        return False

def _make_dendrogram(linkage, names, **kwargs):
    '''
    Currently used by winner plot

    names can be gotten like:
    db = db.pivot("reference","querry","ani")
    names = list(db.columns)
    '''
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
    fig.set_size_inches(10,_x_fig_size(len(names),factor=.2))
    plt.subplots_adjust(left=0.3)

    # Adjust the x labels
    plt.tick_params(axis='both', which='major', labelsize=8)
    axes = plt.gca()
    labels = axes.xaxis.get_majorticklocs()
    for i, label in enumerate(labels):
        labels[i] = (1 - float(label)) * 100
    axes.set_xticklabels(labels)

def _make_special_dendrogram(linkage, names, **kwargs):
    '''
    Make the dendrogram used in plot 2

    names can be gotten like:
        db = db.pivot("reference","querry","ani")
        names = list(db.columns)

    Args:
        linkage: result of scipy.cluster.hierarchy.linkage
        names: names of the linkage

    Kwargs:
        name2cluster: dict
        self_thresh: x-axis for soft line
        threshold: x-axis for hard line
        title_sting: title of the plot
        subtitle_string: subtitle of the plot
        winners: list of "winning" genomes (to be marked with star)
        genome2taxonomy: dictionary to add taxonomy information

    Returns:
        Matplotlib primed with a figure
    '''
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
    fig.set_size_inches(10,_x_fig_size(len(names),factor=.5))
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

def _get_toshow(Gdb):
    '''
    From Gdb, figure out what columns you can show.
    '''
    cols = list(Gdb.columns)
    cols.remove('genome')
    return cols

def _make_scoring_plot(db, bars,**kwargs):
    '''
    Used by winner plot

    db is the database to plot- must contain 'genome' and all columns listed in 'bars'
    bars is all of the columns in the database to become bars
    for taxonomy, put genome2taxonomy in kwargs
    '''
    sns.reset_orig()

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
    fig.set_size_inches(12,_x_fig_size(len(labels), factor=1))
    plt.subplots_adjust(left=0.5)

def _make_mds_plot(name, dist, names, **kwargs):
    '''
    Use MDS to cluster points.

    Based on:
    http://baoilleach.blogspot.com/2014/01/convert-distance-matrix-to-2d.html

    Args:
        name: title of plot
        dist: linkage databases
        names: list of names in linkage database

    Kwargs:
        annotate: if True, write names of all points
        colors: list of colors to use
        shepard: if True, make shepard plot
        tick_spacing: default = .01
        c2c: cluster to color

    Returns:
        Primes plot in matplotlib
    '''

    # load kwargs
    annotate = kwargs.get('annotate', False)
    colors = kwargs.get('colors', False)
    shepard = kwargs.get('shepard', False)
    tick_spacing = kwargs.get('tick_spacing', .01)
    c2c = kwargs.get('cluster2color', False)

    # perform MDS
    mds = manifold.MDS(n_components=2, dissimilarity="precomputed", random_state=6)
    results = mds.fit(dist)
    coords = results.embedding_

    # make shepard plot
    if shepard:
        _shepard_plot(coords, dist, names)

    # plot
    sns.set_style('whitegrid')
    plt.subplots_adjust(bottom = 0.1)
    plt.scatter(
        coords[:, 0], coords[:, 1], marker = 'o', linestyle = 'None',\
        color = colors
        )

    if annotate:
        for label, x, y in zip(names, coords[:, 0], coords[:, 1]):
            plt.annotate(
                label,
                xy = (x, y), xytext = (-20, 20),
                textcoords = 'offset points', ha = 'right', va = 'bottom',
                bbox = dict(boxstyle = 'round,pad=0.5', fc = 'yellow', alpha = 0.5),
                arrowprops = dict(arrowstyle = '->', connectionstyle = 'arc3,rad=0'))

    plt.title('Primary cluster {0}; grid spacing: {1}% ANI'.format(name, tick_spacing*100))

    plt.axis('equal')
    ax = plt.gca()
    ax.set_xticklabels([])
    ax.set_yticklabels([])
    ax.xaxis.set_major_locator(ticker.MultipleLocator(tick_spacing))
    ax.yaxis.set_major_locator(ticker.MultipleLocator(tick_spacing))

    # Add legend
    legendC = c2c
    if legendC != None:
        patches = [mpatches.Patch(color = color, label = label)
                    for label, color in zip(legendC.keys(), legendC.values())]
        lgd = plt.legend(patches, legendC.keys(), loc = 'center left', \
                   bbox_to_anchor=(1, 0.5), prop = {'size': 10, 'style': 'italic'})

def _shepard_plot(coords, dist, names):
    '''
    A componant of the MDS plot
    '''
    table = {'ani_dist':[], 'mds_dist':[]}
    for v1, mx, my in zip(names, coords[:, 0], coords[:, 1]):
        for v2, mx2, my2 in zip(names, coords[:, 0], coords[:, 1]):
            m_dist = calc_dist(mx, my, mx2, my2)

            a_dist = dist.ix[v1, v2]

            table['mds_dist'].append(m_dist)
            table['ani_dist'].append(a_dist)

    adb = pd.DataFrame(table)
    sns.regplot(data=adb, x='ani_dist', y='mds_dist')


'''
****************************    DEPREICATED   *******************************

*   This is where that re-cluster stuff used to be

################################################################################
'''

def cluster_test_wrapper(wd, **kwargs):
    '''
    DEPRICATED
    '''
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
    Xdb = drep.d_cluster.compare_genomes(bdb,comp_method,wd,**kwargs)

    # Remove values without enough coverage
    if comp_method == 'ANIn':
        Xdb.loc[Xdb['alignment_coverage'] <= cov_thresh, 'ani'] = 0

    # Make it symmetrical
    drep.d_cluster.add_avani(Xdb)
    #Xdb['av_ani'] = Xdb.apply(lambda row: drep.d_cluster.average_ani (row,Xdb),axis=1)
    Xdb['dist'] = 1 - Xdb['av_ani']
    db = Xdb.pivot("reference","querry","dist")

    # Cluster it
    if threshold == None:
        threshold = float(0)
    cdb, linkage = drep.d_cluster.cluster_hierarchical(db, linkage_method = clust_method, \
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

def plot_clustertest(linkage, names, wd, **kwargs):
    '''
    DEPREICATED

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
