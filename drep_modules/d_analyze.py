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

# !!! This is just for testing purposes, obviously
import sys
sys.path.append('/home/mattolm/Programs/drep/')
import drep_modules as dm
import drep_modules
import drep_modules.d_cluster as dClust
import drep_modules.d_filter as dFilter

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
    wd = drep_modules.WorkDirectory.WorkDirectory(wd)
        
    if kwargs.get('cluster_visualization') != None:
        logging.info("calling cluster_vis_wrapper")
        cluster_vis_wrapper(wd, **kwargs)
        
    if kwargs.get('mash_cluster') != None:
        logging.info("calling cluster_test_wrapper")
        cluster_test_wrapper(wd, **kwargs)
        
def cluster_test_wrapper(wd, **kwargs):
    
    # Validate arguments
    mash_cluster = kwargs.get('mash_cluster')
    comp_method = kwargs.get('clustering_method','ANIn')
    assert comp_method in ['ANIn','gANI']
    clust_method = kwargs.get('clusterAlg')
    threshold = kwargs.pop('threshold',None)
    if threshold != None: threshold = 1- float(threshold)
    
    # Make a bdb listing the genomes to cluster
    Cdb = wd.get_db('Cdb')
    Bdb = wd.get_db('Bdb')
    genomes = Cdb['genome'][Cdb['MASH_cluster'] == int(mash_cluster)].tolist()
    bdb = Bdb[Bdb['genome'].isin(genomes)]
    
    # Make the db and linkage
    linkage = None
    names = None
    
    # Make the comparison database
    Xdb = compare_genomes(bdb,comp_method,wd,**kwargs)
    
    # Make it symmetrical
    Xdb['av_ani'] = Xdb.apply(lambda row: dClust.average_ani (row,Xdb),axis=1)
    Xdb['dist'] = 1 - Xdb['av_ani']
    db = Xdb.pivot("reference","querry","dist")
    
    # Cluster it
    if threshold == None:
        threshold = float(0)
    cdb, linkage = dClust.cluster_hierarchical(db, linkage_method = clust_method, \
                            linkage_cutoff = threshold)
    
    # Get names
    names = list(db.columns)
    
    # Make the plot directory
    plot_dir = wd.location + '/figures/cluster_tests/'
    if not os.path.exists(plot_dir):
        os.makedirs(plot_dir)
        
    # Make the plot
    plot_single_dendrogram(linkage, names, threshold=threshold, \
        title="{0}_MASH_cluster_{1}".format(comp_method,mash_cluster), loc = plot_dir, \
        **kwargs)
    

def compare_genomes(bdb, comp_method, wd, **kwargs):
    p = kwargs.get('processors')
    genomes = bdb['location'].tolist()
    
    if comp_method == 'ANIn':
        ANIn_folder = wd.location + '/data/ANIn_files/'
        if not os.path.exists(ANIn_folder):
            os.makedirs(ANIn_folder)
            
        # Gen commands
        cmds = []
        files = []
        for g1 in genomes:
            for g2 in genomes:
                file_name = "{0}{1}_vs_{2}".format(ANIn_folder, \
                            dClust.get_genome_name_from_fasta(g1),\
                            dClust.get_genome_name_from_fasta(g2))
                files.append(file_name)
                            
                # If the file doesn't already exist, add it to what needs to be run
                if not os.path.isfile(file_name + '.delta'):
                    cmds.append(dClust.gen_nucmer_cmd(file_name,g1,g2))
                    
        # Run commands
        if len(cmds) > 0:
            dClust.thread_nucmer_cmds_status(cmds,p)
            
        # Parse output
        org_lengths = {dClust.get_genome_name_from_fasta(y):dm.fasta_length(x) for x,y in zip(genomes,genomes)}
        deltafiles = ["{0}.delta".format(file) for file in files]
        df = dClust.process_deltafiles(deltafiles, org_lengths)
        
        return df
        
    elif comp_method == 'gANI':
        gANI_exe = kwargs.get('gANI_exe','/home/mattolm/download/ANIcalculator_v1/ANIcalculator')
        
        # Make folders
        gANI_folder = wd.location + '/data/gANI_files/'
        if not os.path.exists(gANI_folder):
            os.makedirs(gANI_folder)
        # Remove this folder- it shouldn't exist and if it does it messes things up
        crap_folder = gANI_folder + 'ani.blast.dir/'
        if os.path.exists(crap_folder):
            logging.info("CRAP FOLDER EXISTS FOR gANI- removing {0}".format(crap_folder))
            shutil.rmtree(crap_folder)
        
        prod_folder = wd.location + '/data/prodigal/'
        if not os.path.exists(prod_folder):
            os.makedirs(prod_folder)
        
        # Run prodigal
        print("Running prodigal...")
        dFilter.run_prodigal(bdb, prod_folder)
        
        # Gen gANI commands
        print("Running gANI...")
        cmds = []
        files = []
        for i, g1 in enumerate(genomes):
            for j, g2 in enumerate(genomes):
                if i > j:
                    name1= dClust.get_genome_name_from_fasta(g1)
                    name2= dClust.get_genome_name_from_fasta(g2)
                    file_name = "{0}{1}_vs_{2}.gANI".format(gANI_folder, \
                                name1, name2)
                    files.append(file_name)
                            
                    # If the file doesn't already exist, add it to what needs to be run
                    if not os.path.isfile(file_name):
                        fna1 = "{0}{1}.fna".format(prod_folder,name1)
                        fna2 = "{0}{1}.fna".format(prod_folder,name2)
                        cmds.append(dClust.gen_gANI_cmd(file_name,fna1,fna2,gANI_folder,gANI_exe))

                    
        # Run commands
        if len(cmds) > 0:
            logging.info('Running gANI commands: {0}'.format('\n'.join([' '.join(x) for x in cmds])))
            dClust.thread_mash_cmds_status(cmds,p)
        else:
            print("gANI already run- will not re-run")
            
        # Parse output
        df = dClust.process_gani_files(files)
        
        return df

def cluster_vis_wrapper(wd, **kwargs):
    
    # Make the plot directory
    plot_dir = wd.location + '/figures/'
    if not os.path.exists(plot_dir):
        os.makedirs(plot_dir)
        
    # Determine what plots to make
    arg_plots = kwargs.get('cluster_visualization')
    options = ['1','2','3','4','5','6','7']
    to_plot = parse_options(options, arg_plots)
    logging.info("making cviz plots {0}".format(to_plot))
    print("Plotting cviz plots {0}".format(' '.join(sorted(to_plot))))
    
    # 1) MASH cluster-map
    if '1' in to_plot:
         # Load the required data
         Mdb = wd.get_db('Mdb')
         Cdb = wd.get_db('Cdb')
         Mlinkage = wd.get_MASH_linkage()
         
         clust_args = wd.arguments['cluster']
         ML_thresh = clust_args.get('M_Lcutoff', False)
         
         # Make the plot
         print("Plotting MASH clustermap...")
         plot_MASH_clustermap(Mdb, Cdb, Mlinkage, threshold = ML_thresh,\
                            plot_dir = plot_dir)
    
    # 2) MASH dendrogram
    if '2' in to_plot:
        # Load the required data
         Mdb = wd.get_db('Mdb')
         Cdb = wd.get_db('Cdb')
         Mlinkage = wd.get_MASH_linkage()
         
         clust_args = wd.arguments['cluster']
         ML_thresh = clust_args.get('M_Lcutoff', False)
         
         # Make the plot
         print("Plotting MASH dendrogram...")
         plot_MASH_dendrogram(Mdb, Cdb, Mlinkage, threshold = ML_thresh,\
                            plot_dir = plot_dir)
                            
    # 3) ANIn clustermaps
    if '3' in to_plot:
        # Load the required data
         Ndb = wd.get_db('Ndb')
         Cdb = wd.get_db('Cdb')
         Nlinkage = wd.get_ANIn_linkages()
         
         clust_args = wd.arguments['cluster']
         NL_thresh = clust_args.get('N_Lcutoff', False)
         
         # Make the plot
         print("Plotting ANIn clustermaps...")
         plot_ANIn_clustermaps(Ndb, Cdb, Nlinkage, threshold = NL_thresh,\
                            plot_dir = plot_dir)
                            
    # 4) ANIn dendrograms
    if '4' in to_plot:
        # Load the required data
         Ndb = wd.get_db('Ndb')
         Cdb = wd.get_db('Cdb')
         Nlinkage = wd.get_ANIn_linkages()
         
         clust_args = wd.arguments['cluster']
         NL_thresh = clust_args.get('N_Lcutoff', False)
         
         # Make the plot
         print("Plotting ANIn dendrograms...")
         plot_ANIn_dendrograms(Ndb, Cdb, Nlinkage, threshold = NL_thresh,\
                            plot_dir = plot_dir)
                            
   # 5) Comparison scatterplots
    if '5' in to_plot:
        # Load the required data
         Ndb = wd.get_db('Ndb')
         Mdb = wd.get_db('Mdb')
         
         # Make the plot
         print("Plotting Scatterplots...")
         plot_scatterplots(Mdb, Ndb, plot_dir = plot_dir)
         
    # 6) Simple bin scorring
    if '6' in to_plot:
        # Load the required data
        Sdb = wd.get_db('Sdb')
        Cdb = wd.get_db('Cdb')
        Wdb = wd.get_db('Wdb')
            
        # Make the plot
        print("Plotting simple scorring plot...")
        plot_winner_scoring(Wdb, Sdb, Cdb, plot_dir = plot_dir)
        
    # 7) Complex bin scorring
    if '7' in to_plot:
        # Load the required data
        Sdb = wd.get_db('Sdb')
        Cdb = wd.get_db('Cdb')
        Wdb = wd.get_db('Wdb')
        Chdb = wd.get_db('Chdb')
            
        # Make the plot
        print("Plotting complex scorring plot...")
        plot_winner_scoring_complex(Wdb, Sdb, Cdb, Chdb, plot_dir = plot_dir)

    
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
        pp = PdfPages(plot_dir + 'Cviz_scatterplots.pdf')
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
    
def plot_MASH_dendrogram(Mdb, Cdb, linkage, threshold = False, plot_dir = False, **kwargs):    
    db = Mdb.pivot("genome1","genome2","similarity")
    names = list(db.columns)
    name2cluster = Cdb.set_index('genome')['MASH_cluster'].to_dict()
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
        plt.savefig(plot_dir + 'Cviz_MASH_dendrogram.pdf')
    plt.show()
    plt.close('all')
    
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
        name2cluster = Cdb.set_index('genome')['ANIn_cluster'].to_dict()
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
        c_genomes = Cdb['genome'][Cdb['MASH_cluster'] == int(cluster)]
        db = Ndb[Ndb['reference'].isin(c_genomes)]
        
        # Get the highest self-comparison
        self_thresh = 1 - db['ani'][db['reference'] == db['querry']].max()
        if self_thresh == float(0):  self_thresh = 1.0e-7 # Still draw if 0
       
        # Get the colors set up
        db = db.pivot("reference","querry","ani")
        names = list(db.columns)
        name2cluster = Cdb.set_index('genome')['ANIn_cluster'].to_dict()
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
def plot_winner_scoring(Wdb, Sdb, Cdb, plot_dir= False, **kwargs):
    save = False
    
    # Initialize a .pdf
    if plot_dir != False:
        pp = PdfPages(plot_dir + 'Cviz_cluster_scoring.pdf')
        save = True
        
    winners = list(Wdb['genome'].unique())
    for cluster in sorted(Cdb['ANIn_cluster'].unique(), key=lambda x: comp_cluster(x)):
        d = Cdb[Cdb['ANIn_cluster'] == cluster]
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
        pp = PdfPages(plot_dir + 'Cviz_cluster_scoring_complex.pdf')
        save = True
    
    # Figure out what you're going to show    
    winners = list(Wdb['genome'].unique())
    Chdb['genome'] = Chdb['Bin Id']
    bars = kwargs.get('to_show',['score','N50 (scaffolds)','Completeness','Contamination','Genome size (bp)', 'Strain heterogeneity'])
    
    for cluster in sorted(Cdb['ANIn_cluster'].unique(), key=lambda x: comp_cluster(x)):
        # Make a db for this cluster
        d = Cdb[Cdb['ANIn_cluster'] == cluster]
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
            g.annotate("{0:.1f}".format(vals[i]), (p.get_width(), p.get_y()+(p.get_height()/2) ))
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
        
        fig = plt.gcf()
        fig.set_size_inches(10,x_fig_size(len(labels), factor=1))
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
        max_value = df[feature_name].max()
        result[feature_name] = df[feature_name] / max_value 
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
    cm = plt.get_cmap('gist_rainbow')
    cm = rand_cmap(len(set(name2cluster.values())))
    
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