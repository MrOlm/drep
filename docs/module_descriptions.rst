Module Descriptions
===================

The functionality of dRep is broken up into modules. The user can run the modules separately, or together in workflows. For example, you could run::

 $ dRep filter example_workD -g path/to/genomes*.fasta

 $ dRep cluster example_workD

 $ dRep analyze example_workD -pl a

OR::

 $ dRep comparison_wf example_workD -g path/to/genomes*.fasta

There are two ways of doing the same thing. To see a list of available modules, check the help::

  mattolm@Matts-MacBook-Pro:~/Programs/drep/docs$ dRep -h

                  ...::: dRep v0.2.0 :::...

    Choose one of the operations below for more detailed help.
    Example: dRep dereplicate_wf -h

    Workflows:
      dereplicate_wf  -> Combine several of the operations below to de-replicate a genome list
      compare_wf      -> Simply compare a list of genomes

    Single opterations:
      filter          -> Filter a genome list based on size, completeness, and/or contamination
      cluster         -> Compare and cluster a genome list based on MASH and ANIn/gANI
      adjust          -> Adjust genome clusters
      choose          -> Choose the best genome from each genome cluster
      evaluate        -> Evaluate genome de-replication
      bonus           -> Other random operations (currently just determine taxonomy)
      analyze         -> Make figures realted to the above operations; test alternative clustering

Work Directory
--------------

The work directory is where all of the program's internal workings, log files, cached data, and output is stored. When running dRep modules multiple times on the same dataset, **it is essential** that you use the same work directory so the program can find the results of previous runs.

.. seealso::

  :doc:`interpreting_output`
    for help finding where the output from your run is located in the work directory

  :doc:`advanced_use`
    for access to the raw internal data (which can be very useful)

Filter
------

Filter is used filter the genome set (for why this is necessary, see :doc:`choosing_parameters`). This is done using checkM. All genomes which don't pass the length threshold are filtered first to avoid running checkM unnecessarily. All genomes which don't pass checkM thresholds are filtered before comparisons are run to avoid running comparisons unnecessarily.

.. warning::

  Due to a bug in checkM, all genomes must have at least one ORF called or else checkM will stall. So a length minimum of at least 10,000bp is recommended.

To see the command-line options, check the help::

  mattolm@Matts-MacBook-Pro:~/Programs/drep/docs$ dRep filter -h
  usage: dRep filter [-p PROCESSORS] [-d] [-o] [-h] [-l LENGTH]
                     [-comp COMPLETENESS] [-con CONTAMINATION] [-str STRAIN_HTR]
                     [--skipCheckM] [-g [GENOMES [GENOMES ...]]] [--Chdb CHDB]
                     [--checkM_method {lineage_wf,taxonomy_wf}]
                     work_directory

  positional arguments:
    work_directory        Directory where data and output
                          *** USE THE SAME WORK DIRECTORY FOR ALL DREP OPERATIONS ***

  SYSTEM PARAMETERS:
    -p PROCESSORS, --processors PROCESSORS
                          threads (default: 6)
    -d, --dry             dry run- dont do anything (default: False)
    -o, --overwrite       overwrite existing data in work folder (default:
                          False)
    -h, --help            show this help message and exit

  FILTERING OPTIONS:
    -l LENGTH, --length LENGTH
                          Minimum genome length (default: 500000)
    -comp COMPLETENESS, --completeness COMPLETENESS
                          Minumum genome completeness (default: 75)
    -con CONTAMINATION, --contamination CONTAMINATION
                          Maximum genome contamination (default: 25)
    -str STRAIN_HTR, --strain_htr STRAIN_HTR
                          Maximum strain heterogeneity (default: 25)
    --skipCheckM          Don't run checkM- will ignore con and comp settings
                          (default: False)

  I/O PARAMETERS:
    -g [GENOMES [GENOMES ...]], --genomes [GENOMES [GENOMES ...]]
                          genomes to filter in .fasta format. Not necessary if
                          Bdb or Wdb already exist (default: None)
    --Chdb CHDB           checkM run already completed. Must be in --tab_table
                          format. (default: None)
    --checkM_method {lineage_wf,taxonomy_wf}
                          Either lineage_wf (more accurate) or taxonomy_wf
                          (faster) (default: lineage_wf)

Cluster
-------

Cluster is the module that does the actual primary and secondary comparisons. Choosing parameters here can get a bit complicated- see :doc:`choosing_parameters` for information.

To see the command-line options, check the help::

  mattolm@Matts-MacBook-Pro:~/Programs/drep/docs$ dRep cluster -h
  usage: dRep cluster [-p PROCESSORS] [-d] [-o] [-h] [-ms MASH_SKETCH]
                      [-pa P_ANI] [--S_algorithm {ANIn,gANI}] [-sa S_ANI]
                      [-nc COV_THRESH] [-n_PRESET {normal,tight}]
                      [--clusterAlg CLUSTERALG] [--SkipMash] [--SkipSecondary]
                      [-g [GENOMES [GENOMES ...]]]
                      work_directory

  positional arguments:
    work_directory        Directory where data and output
                          *** USE THE SAME WORK DIRECTORY FOR ALL DREP OPERATIONS ***

  SYSTEM PARAMETERS:
    -p PROCESSORS, --processors PROCESSORS
                          threads (default: 6)
    -d, --dry             dry run- dont do anything (default: False)
    -o, --overwrite       overwrite existing data in work folder (default:
                          False)
    -h, --help            show this help message and exit

  CLUSTERING PARAMETERS:
    -ms MASH_SKETCH, --MASH_sketch MASH_SKETCH
                          MASH sketch size (default: 1000)
    -pa P_ANI, --P_ani P_ANI
                          ANI threshold to form primary (MASH) clusters
                          (default: 0.9)
    --S_algorithm {ANIn,gANI}
                          Algorithm for secondary clustering comaprisons
                          (default: ANIn)
    -sa S_ANI, --S_ani S_ANI
                          ANI threshold to form secondary clusters (default:
                          0.99)
    -nc COV_THRESH, --cov_thresh COV_THRESH
                          Minmum level of overlap between genomes when doing
                          secondary comparisons (default: 0.1)
    -n_PRESET {normal,tight}
                          Presents to pass to nucmer
                          tight   = only align highly conserved regions
                          normal  = default ANIn parameters (default: normal)
    --clusterAlg CLUSTERALG
                          Algorithm used to cluster genomes (passed to
                          scipy.cluster.hierarchy.linkage (default: average)
    --SkipMash            Skip MASH clustering, just do secondary clustering on
                          all genomes (default: False)
    --SkipSecondary       Skip secondary clustering, just perform MASH
                          clustering (default: False)

  I/O PARAMETERS:
    -g [GENOMES [GENOMES ...]], --genomes [GENOMES [GENOMES ...]]
                          genomes to cluster in .fasta format. Not necessary if
                          already loaded sequences with the "filter" operation
                          (default: None)

Choose
------

Choose is the module that picks the best genome from each secondary cluster identified in **Cluster**. It does this based off of the formula:

.. math:: score = A(completeness) + B(log_{10}(N50)) – C(contamination) – D(strain heterogeneity) + E(log_{10}(genome size))

Where A-E are command-line arguments, and the genome with the highest score is the "best". By default, A-E are 1,0.5,5,1,0, respectively.

To see the command-line options, check the help::

  mattolm@Matts-MacBook-Pro:~/Programs/drep/docs$ dRep choose -h
  usage: dRep choose [-p PROCESSORS] [-d] [-o] [-h] [-comW COMPLETENESS_WEIGHT]
                     [-conW CONTAMINATION_WEIGHT] [-N50W N50_WEIGHT]
                     [-sizeW SIZE_WEIGHT] [-strW STRAIN_HETEROGENEITY_WEIGHT]
                     [-c CLUSTER] [-t THRESHOLD] [-m {gANI,ANIn}]
                     [-mc MINIMUM_COVERAGE]
                     [-a {average,complete,weighted,single}]
                     [--checkM_method {taxonomy_wf,lineage_wf}]
                     work_directory

  positional arguments:
    work_directory        Directory where data and output
                          *** USE THE SAME WORK DIRECTORY FOR ALL DREP OPERATIONS ***

  SYSTEM PARAMETERS:
    -p PROCESSORS, --processors PROCESSORS
                          threads (default: 6)
    -d, --dry             dry run- dont do anything (default: False)
    -o, --overwrite       overwrite existing data in work folder (default:
                          False)
    -h, --help            show this help message and exit

  SCORRING CHRITERIA
  Based off of the formula: Completeness - Contamination + log(N50) + log(size):
    -comW COMPLETENESS_WEIGHT, --completeness_weight COMPLETENESS_WEIGHT
                          completeness weight (default: 1)
    -conW CONTAMINATION_WEIGHT, --contamination_weight CONTAMINATION_WEIGHT
                          contamination weight (default: 5)
    -N50W N50_WEIGHT, --N50_weight N50_WEIGHT
                          weight of log(genome N50) (default: 0.5)
    -sizeW SIZE_WEIGHT, --size_weight SIZE_WEIGHT
                          weight of log(genome size) (default: 0)
    -strW STRAIN_HETEROGENEITY_WEIGHT, --strain_heterogeneity_weight STRAIN_HETEROGENEITY_WEIGHT
                          strain heterogeneity weight (default: 1)

  RE-CLUSTER PRIMARY CLUSETERS:
    -c CLUSTER, --cluster CLUSTER
                          primary cluster to be adjusted (default: None)
    -t THRESHOLD, --threshold THRESHOLD
                          clustering threshold to apply (default: 0.99)
    -m {gANI,ANIn}, --clustering_method {gANI,ANIn}
                          Clustering method to apply (default: ANIn)
    -mc MINIMUM_COVERAGE, --minimum_coverage MINIMUM_COVERAGE
                          Minimum coverage for ANIn (default: 0.1)
    -a {average,complete,weighted,single}, --clusterAlg {average,complete,weighted,single}
                          Algorithm used to cluster genomes (passed to
                          scipy.cluster.hierarchy.linkage) (default: average)

  OTHER:
    --checkM_method {taxonomy_wf,lineage_wf}
                          Either lineage_wf (more accurate) or taxonomy_wf
                          (faster) (default: lineage_wf)

Analyze
-------

Analyze is the module that makes all of the figures. It also has the option to visualize how a secondary cluster would look with different parameters (for example, using ANIm instead of gANI). To do that, use the ``RE-CLUSTER PRIMARY CLUSTERS`` arguments. To make plots, just use the -pl argument.

To see the command-line options, check the help::

  mattolm@Matts-MacBook-Pro:~/Programs/drep/docs$ dRep analyze -h
  usage: dRep analyze [-p PROCESSORS] [-d] [-o] [-h] [-c CLUSTER] [-t THRESHOLD]
                      [-m {ANIn,gANI}] [-mc MINIMUM_COVERAGE]
                      [-a {complete,average,single,weighted}]
                      [-pl [PLOTS [PLOTS ...]]]
                      work_directory

  positional arguments:
    work_directory        Directory where data and output
                          *** USE THE SAME WORK DIRECTORY FOR ALL DREP OPERATIONS ***

  SYSTEM PARAMETERS:
    -p PROCESSORS, --processors PROCESSORS
                          threads (default: 6)
    -d, --dry             dry run- dont do anything (default: False)
    -o, --overwrite       overwrite existing data in work folder (default:
                          False)
    -h, --help            show this help message and exit

  RE-CLUSTER PRIMARY CLUSETERS:
    -c CLUSTER, --cluster CLUSTER
                          primary cluster to be adjusted (default: None)
    -t THRESHOLD, --threshold THRESHOLD
                          clustering threshold to apply (default: 0.99)
    -m {ANIn,gANI}, --clustering_method {ANIn,gANI}
                          Clustering method to apply (default: ANIn)
    -mc MINIMUM_COVERAGE, --minimum_coverage MINIMUM_COVERAGE
                          Minimum coverage for ANIn (default: 0.1)
    -a {complete,average,single,weighted}, --clusterAlg {complete,average,single,weighted}
                          Algorithm used to cluster genomes (passed to
                          scipy.cluster.hierarchy.linkage) (default: average)

  PLOTTING:
    -pl [PLOTS [PLOTS ...]], --plots [PLOTS [PLOTS ...]]
                          Plots. Input 'all' or 'a' to plot all
                          1) Primary clustering dendrogram
                          2) Secondary clustering dendrograms
                          3) Secondary clusters heatmaps
                          4) Comparison scatterplots
                          5) Cluster scorring plot
                          6) Winning genomes
                           (default: None)

Evaluate
--------

Evaluate performs a series of checks to alert the user to potential problems with de-replication. It has two things that it can look for:

**de-replicated genome similarity**- this is comparing all of the de-replicated genomes to each other and making sure they're not too similar. This is to try and catch cases where similar genomes were split into different primary clusters, and thus failed to be de-replicated. *Depending on the number of de-replicated genomes, this can take a while*

**secondary clusters that were almost different**- this alerts you to cases where genomes are on the edge between being considered "same" or "different", depending on the clustering parameters you used. *This module reads the parameters you used during clustering from the work directory, so you don't need to specify them again.*

To see the command-line options, check the help::

  mattolm@Matts-MacBook-Pro:~/Programs/drep/docs$ dRep evaluate -h
  usage: dRep evaluate [-p PROCESSORS] [-d] [-o] [-h] [--warn_dist WARN_DIST]
                       [--warn_sim WARN_SIM] [--warn_aln WARN_ALN]
                       [-e [EVALUATE [EVALUATE ...]]]
                       work_directory

  positional arguments:
    work_directory        Directory where data and output
                          *** USE THE SAME WORK DIRECTORY FOR ALL DREP OPERATIONS ***

  SYSTEM PARAMETERS:
    -p PROCESSORS, --processors PROCESSORS
                          threads (default: 6)
    -d, --dry             dry run- dont do anything (default: False)
    -o, --overwrite       overwrite existing data in work folder (default:
                          False)
    -h, --help            show this help message and exit

  WARNINGS:
    --warn_dist WARN_DIST
                          How far from the threshold to throw cluster warnings
                          (default: 0.25)
    --warn_sim WARN_SIM   Similarity threshold for warnings between dereplicated
                          genomes (default: 0.98)
    --warn_aln WARN_ALN   Minimum aligned fraction for warnings between
                          dereplicated genomes (ANIn) (default: 0.25)

  EVALUATIONS:
    -e [EVALUATE [EVALUATE ...]], --evaluate [EVALUATE [EVALUATE ...]]
                          Things to evaluate Input 'all' or 'a' to evaluate all
                          1) Evaluate de-replicated genome similarity
                          2) Throw warnings for clusters that were almost different
                          3) Generate a database of information on winning genomes
                           (default: None)

Other
-----

The other modules, **adjust** and **bonus**, are not part of the normal de-replication pipeline but can be very useful.

**adjust** allows the user to change the secondary clustering settings for a single primary cluster. This can be especially helpful when following up on a warning (generated using **evaluate**) to change the way the cluster is made.

**bonus** consists of operations that don't really fit in with the functions of dRep, but can be helpful. Currently the only thing it can do is determine taxonomy of your bins. This is done using centrifuge, similar to how `anvi'o does it <http://merenlab.org/2016/06/18/importing-taxonomy/>`_. If you choose to use this option, the taxonomy of genome will be shown with the filename in most figures.
