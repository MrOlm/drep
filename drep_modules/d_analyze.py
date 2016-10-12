#!/usr/bin/env python3

import matplotlib
matplotlib.use('Agg')

import pandas as pd
import os
import seaborn as sns
from matplotlib import pyplot as plt
import scipy.cluster.hierarchy

"""#############################################################################
                            MODULE ARCHITECTURE

****************************    Clustering   *******************************
Graphs - Heatmaps

*   plot_Mdb_heatmap(Mdb)
    -   Make a single heat-map of Mdb
    
*   plot_ANIn_heatmap(Ndb)
    -   Make a heat-map of ANIn for every MASH cluster
    
*   plot_ANIn_cov_heatmap(Ndb)
    -   Make a heat-map of ANI_alignment_coverage for every MASH cluster
    
*   plot_ANIn_clusters(Ndb, Cdb)
    -   Make a heat-map for every Ndb cluster
    
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
    
Graphs - Custom

*   plot_cluster_tightness(Ndb)
    - Come up with some way of visualizing the variation within ANIn clusters,
      versus the variation between clusters
      - Show the average and max tightness within and between all clusters
################################################################################
"""

METHOD = 'single'


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

def plot_MASH_vs_ANIn_ani(Mdb,Ndb,exclude_zero_MASH=True):
    mdb = Mdb.copy()
    mdb.rename(columns={'genome1':'querry','genome2':'reference',
                        'similarity':'MASH_ANI'},inplace=True)
    if exclude_zero_MASH:
        mdb= mdb[mdb['MASH_ANI'] > 0]
        
    db = pd.merge(mdb,Ndb)
    db.rename(columns={'ani':'ANIn'},inplace=True)
    g = sns.jointplot(x='ANIn',y='MASH_ANI',data=db)
    return g

def plot_MASH_vs_ANIn_cov(Mdb,Ndb,exclude_zero_MASH=True):
    mdb = Mdb.copy()
    mdb.rename(columns={'genome1':'querry','genome2':'reference',
                        'similarity':'MASH_ANI'},inplace=True)
    if exclude_zero_MASH:
        mdb= mdb[mdb['MASH_ANI'] > 0]
        
    db = pd.merge(mdb,Ndb)
    db.rename(columns={'alignment_coverage':'ANIn_cov'},inplace=True)
    g = sns.jointplot(x='ANIn_cov',y='MASH_ANI',data=db)
    return g
    
def plot_ANIn_vs_ANIn_cov(Ndb):
    db = Ndb.copy()
    db.rename(columns={'alignment_coverage':'ANIn_cov','ani':'ANIn'},inplace=True)
    g = sns.jointplot(x='ANIn_cov',y='ANIn',data=db)
    return g
    
"""
CLUSETER PLOTS
"""

def plot_MASH_clusters(Mdb, Cdb, linkage, threshold= False):
    gs = list()
    
    db = Mdb.pivot("genome1","genome2","similarity")
    names = list(db.columns)
    name2cluster = Cdb.set_index('genome')['MASH_cluster'].to_dict()
    colors = gen_color_list(names, name2cluster)
    name2color = gen_color_dictionary(names, name2cluster)
    
    
    # Make the dendrogram
    g = fancy_dengrogram(linkage,names,name2color,threshold=0.1)
    plt.title('MASH clustering')
    plt.ylabel('distance (about 1 - MASH_ANI)')
    fig = plt.gcf()
    gs.append(fig)
    plt.clf()
    
    
    # Make the clustermap
    g = sns.clustermap(db, row_linkage = linkage, col_linkage = linkage, \
                        row_colors = colors, col_colors = colors)
    g.fig.suptitle("MASH clustering")
    plt.setp(g.ax_heatmap.yaxis.get_majorticklabels(), rotation=0)
    fig = plt.gcf()
    gs.append(fig)
    plt.clf()
    
    return gs
    

def plot_ANIn_clusters(Ndb, Cdb, cluster2linkage, threshold= False):
    gs = list()
    
    for cluster in cluster2linkage.keys():
        linkage = cluster2linkage[cluster]
        
        # Refine Ndb to just have the clusters of the linkage
        c_genomes = Cdb['genome'][Cdb['MASH_cluster'] == int(cluster)]
        db = Ndb[Ndb['reference'].isin(c_genomes)]
        db = db.pivot("reference","querry","ani")
        
        # Get the colors set up
        names = list(db.columns)
        name2cluster = Cdb.set_index('genome')['ANIn_cluster'].to_dict()
        colors = gen_color_list(names, name2cluster)
        name2color = gen_color_dictionary(names, name2cluster)
    
        # Make the dendrogram
        g = fancy_dengrogram(linkage,names,name2color,threshold=0.1)
        plt.title('ANI of MASH cluster {0}'.format(cluster))
        plt.ylabel('distance (about 1 - ANIn)')
        fig = plt.gcf()
        plt.close()
        gs.append(fig)
        gs.append(g)
    
        # Make the clustermap
        g = sns.clustermap(db, row_linkage = linkage, col_linkage = linkage, \
                            row_colors = colors, col_colors = colors)
        g.fig.suptitle('ANI of MASH cluster {0}'.format(cluster))
        plt.setp(g.ax_heatmap.yaxis.get_majorticklabels(), rotation=0)
        fig = plt.gcf()
        plt.close()
        gs.append(fig)
    
    return gs
    
"""
OTHER
"""

def fancy_dengrogram(linkage,names,name2color,threshold=False):
    
    # Make the dendrogram
    scipy.cluster.hierarchy.dendrogram(linkage,labels=names)
    plt.xticks(rotation=90)
    
    # Color the names
    ax = plt.gca()
    xlbls = ax.get_xmajorticklabels()
    for lbl in xlbls:
        lbl.set_color(name2color[lbl.get_text()])
        
    # Add the threshold
    if threshold:
        plt.axhline(y=threshold, c='k')
        
    g = plt.gcf()
    return g

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

def test_clustering():
    print("You should make some test cases here!")

if __name__ == '__main__':
	test_clustering()