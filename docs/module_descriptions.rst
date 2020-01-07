Module Descriptions
===================

The functionality of dRep is broken up into modules. The user can run the modules separately, or together in workflows. For example, you could run::

 $ dRep filter example_workD -g path/to/genomes*.fasta

 $ dRep cluster example_workD

 $ dRep analyze example_workD -pl a

OR::

 $ dRep compare example_workD -g path/to/genomes*.fasta

There are two ways of doing the same thing. To see a list of available modules, check the help::

 $ dRep -h

                ...::: dRep v2.0.0 :::...

  Choose one of the operations below for more detailed help.
  Example: dRep dereplicate -h

  Workflows:
    dereplicate  -> Combine several of the operations below to de-replicate a genome list
    compare      -> Simply compare a list of genomes

  Single operations:
    filter          -> Filter a genome list based on size, completeness, and/or contamination
    cluster         -> Compare and cluster a genome list based on MASH and ANIn/gANI
    choose          -> Choose the best genome from each genome cluster
    evaluate        -> Evaluate genome de-replication
    bonus           -> Other random operations (currently just determine taxonomy)
    analyze         -> Make figures related to the above operations; test alternative clustering

Work Directory
--------------

The work directory is where all of the program's internal workings, log files, cached data, and output is stored. When running dRep modules multiple times on the same dataset, **it is essential** that you use the same work directory so the program can find the results of previous runs.

.. seealso::

  :doc:`example_output`
    for help finding where the output from your run is located in the work directory

  :doc:`advanced_use`
    for access to the raw internal data (which can be very useful)

Compare and Dereplicate
------
These are higher-level operations that call the modules below in succession.

Compare runs the modules:

* cluster
* bonus
* evaluate
* analyze

Dereplicate runs the modules:

* filter
* cluster
* choose
* bonus
* evaluate
* analyze

Filter
------

Filter is used filter the genome set (for why this is necessary, see :doc:`choosing_parameters`). This is done using checkM. All genomes which don't pass the length threshold are filtered first to avoid running checkM unnecessarily. All genomes which don't pass checkM thresholds are filtered before comparisons are run to avoid running comparisons unnecessarily.

.. warning::

  All genomes must have at least one ORF called or else checkM will stall, so a length minimum of at least 10,000bp is recommended.

To see the command-line options, check the help::

  $ dRep filter -h
  usage: dRep filter [-p PROCESSORS] [-d] [-h] [-l LENGTH] [-comp COMPLETENESS]
                   [-con CONTAMINATION] [--ignoreGenomeQuality]
                   [-g [GENOMES [GENOMES ...]]] [--genomeInfo GENOMEINFO]
                   [--checkM_method {taxonomy_wf,lineage_wf}]
                   [--set_recursion SET_RECURSION]
                   work_directory

  positional arguments:
  work_directory        Directory where data and output
                        *** USE THE SAME WORK DIRECTORY FOR ALL DREP OPERATIONS ***

  SYSTEM PARAMETERS:
  -p PROCESSORS, --processors PROCESSORS
                        threads (default: 6)
  -d, --debug           make extra debugging output (default: False)
  -h, --help            show this help message and exit

  FILTERING OPTIONS:
  -l LENGTH, --length LENGTH
                        Minimum genome length (default: 50000)
  -comp COMPLETENESS, --completeness COMPLETENESS
                        Minumum genome completeness (default: 75)
  -con CONTAMINATION, --contamination CONTAMINATION
                        Maximum genome contamination (default: 25)
  --ignoreGenomeQuality
                        Don't run checkM or do any quality filtering. NOT
                        RECOMMENDED! This is useful for use with
                        bacteriophages or eukaryotes or things where checkM
                        scoring does not work. Will only choose genomes based
                        on length and N50 (default: False)

  I/O PARAMETERS:
  -g [GENOMES [GENOMES ...]], --genomes [GENOMES [GENOMES ...]]
                        genomes to filter in .fasta format. Not necessary if
                        Bdb or Wdb already exist (default: None)
  --genomeInfo GENOMEINFO
                        location of .csv file containing quality information
                        on the genomes. Must contain: ["genome"(basename of
                        .fasta file of that genome), "completeness"(0-100
                        value for completeness of the genome),
                        "contamination"(0-100 value of the contamination of
                        the genome)] (default: None)
  --checkM_method {taxonomy_wf,lineage_wf}
                        Either lineage_wf (more accurate) or taxonomy_wf
                        (faster) (default: lineage_wf)
  --set_recursion SET_RECURSION
                        Increases the python recursion limit. NOT RECOMMENDED
                        unless checkM is crashing due to recursion issues.
                        Recommended to set to 2000 if needed, but setting this
                        could crash python (default: 0)

Cluster
-------

Cluster is the module that does the actual primary and secondary comparisons. Choosing parameters here can get a bit complicated- see :doc:`choosing_parameters` for information.

To see the command-line options, check the help::

  $ dRep cluster -h
  usage: dRep cluster [-p PROCESSORS] [-d] [-h] [-ms MASH_SKETCH]
                    [--S_algorithm {ANIn,gANI,ANImf,goANI}]
                    [-n_PRESET {normal,tight}] [-pa P_ANI] [-sa S_ANI]
                    [--SkipMash] [--SkipSecondary] [-nc COV_THRESH]
                    [-cm {total,larger}] [--clusterAlg CLUSTERALG]
                    [-g [GENOMES [GENOMES ...]]]
                    work_directory

  positional arguments:
  work_directory        Directory where data and output
                        *** USE THE SAME WORK DIRECTORY FOR ALL DREP OPERATIONS ***

  SYSTEM PARAMETERS:
  -p PROCESSORS, --processors PROCESSORS
                        threads (default: 6)
  -d, --debug           make extra debugging output (default: False)
  -h, --help            show this help message and exit

  GENOME COMPARISON PARAMETERS:
  -ms MASH_SKETCH, --MASH_sketch MASH_SKETCH
                        MASH sketch size (default: 1000)
  --S_algorithm {ANIn,gANI,ANImf,goANI}
                        Algorithm for secondary clustering comaprisons:
                        ANImf = (RECOMMENDED) Align whole genomes with nucmer; filter alignment; compare aligned regions
                        ANIn  = Align whole genomes with nucmer; compare aligned regions
                        gANI  = Identify and align ORFs; compare aligned ORFS
                         (default: ANImf)
  -n_PRESET {normal,tight}
                        Presets to pass to nucmer
                        tight   = only align highly conserved regions
                        normal  = default ANIn parameters (default: normal)

  CLUSTERING PARAMETERS:
  -pa P_ANI, --P_ani P_ANI
                        ANI threshold to form primary (MASH) clusters
                        (default: 0.9)
  -sa S_ANI, --S_ani S_ANI
                        ANI threshold to form secondary clusters (default:
                        0.99)
  --SkipMash            Skip MASH clustering, just do secondary clustering on
                        all genomes (default: False)
  --SkipSecondary       Skip secondary clustering, just perform MASH
                        clustering (default: False)
  -nc COV_THRESH, --cov_thresh COV_THRESH
                        Minmum level of overlap between genomes when doing
                        secondary comparisons (default: 0.1)
  -cm {total,larger}, --coverage_method {total,larger}
                        Method to calculate coverage of an alignment
                        (for ANIn/ANImf only; gANI can only do larger method)
                        total   = 2*(aligned length) / (sum of total genome lengths)
                        larger  = max((aligned length / genome 1), (aligned_length / genome2))
                         (default: larger)
  --clusterAlg CLUSTERALG
                        Algorithm used to cluster genomes (passed to
                        scipy.cluster.hierarchy.linkage (default: average)

  I/O PARAMETERS:
  -g [GENOMES [GENOMES ...]], --genomes [GENOMES [GENOMES ...]]
                        genomes to cluster in .fasta format. Not necessary if
                        already loaded sequences with the "filter" operation
                        (default: None)

Choose
------

Choose is the module that picks the best genome from each secondary cluster identified in **Cluster**. It does this based off of the formula:

.. math:: score = A(completeness) – B(contamination) +  C(Contamination * (strain_heterogeneity/100)) + D(log(N50)) + E(log(size))

Where A-E are command-line arguments, and the genome with the highest score is the "best". By default, A-E are 1,5,1,0.5,0 respectively.

To see the command-line options, check the help::

  $ dRep choose -h
  usage: dRep choose [-p PROCESSORS] [-d] [-h] [-comW COMPLETENESS_WEIGHT]
                   [-conW CONTAMINATION_WEIGHT]
                   [-strW STRAIN_HETEROGENEITY_WEIGHT] [-N50W N50_WEIGHT]
                   [-sizeW SIZE_WEIGHT]
                   [--checkM_method {lineage_wf,taxonomy_wf}]
                   [--genomeInfo GENOMEINFO] [--ignoreGenomeQuality]
                   work_directory

  positional arguments:
  work_directory        Directory where data and output
                        *** USE THE SAME WORK DIRECTORY FOR ALL DREP OPERATIONS ***

  SYSTEM PARAMETERS:
  -p PROCESSORS, --processors PROCESSORS
                        threads (default: 6)
  -d, --debug           make extra debugging output (default: False)
  -h, --help            show this help message and exit

  SCORING CRITERIA
  Based off of the formula:
  A*Completeness - B*Contamination + C*(Contamination * (strain_heterogeneity/100)) + D*log(N50) + E*log(size)

  A = completeness_weight; B = contamination_weight; C = strain_heterogeneity_weight; D = N50_weight; E = size_weight:
  -comW COMPLETENESS_WEIGHT, --completeness_weight COMPLETENESS_WEIGHT
                        completeness weight (default: 1)
  -conW CONTAMINATION_WEIGHT, --contamination_weight CONTAMINATION_WEIGHT
                        contamination weight (default: 5)
  -strW STRAIN_HETEROGENEITY_WEIGHT, --strain_heterogeneity_weight STRAIN_HETEROGENEITY_WEIGHT
                        strain heterogeneity weight (default: 1)
  -N50W N50_WEIGHT, --N50_weight N50_WEIGHT
                        weight of log(genome N50) (default: 0.5)
  -sizeW SIZE_WEIGHT, --size_weight SIZE_WEIGHT
                        weight of log(genome size) (default: 0)

  OTHER:
  --checkM_method {lineage_wf,taxonomy_wf}
                        Either lineage_wf (more accurate) or taxonomy_wf
                        (faster) (default: lineage_wf)
  --genomeInfo GENOMEINFO
                        location of .csv file containing quality information
                        on the genomes. Must contain: ["genome"(basename of
                        .fasta file of that genome), "completeness"(0-100
                        value for completeness of the genome),
                        "contamination"(0-100 value of the contamination of
                        the genome)] (default: None)
  --ignoreGenomeQuality
                        Don't run checkM or do any quality filtering. NOT
                        RECOMMENDED! This is useful for use with
                        bacteriophages or eukaryotes or things where checkM
                        scoring does not work. Will only choose genomes based
                        on length and N50 (default: False)

Analyze
-------

Analyze is the module that makes all of the figures.

To see the command-line options, check the help::

  $ dRep analyze -h
  usage: dRep analyze [-p PROCESSORS] [-d] [-h] [-pl [PLOTS [PLOTS ...]]]
                      work_directory

  positional arguments:
    work_directory        Directory where data and output
                          *** USE THE SAME WORK DIRECTORY FOR ALL DREP OPERATIONS ***

  SYSTEM PARAMETERS:
    -p PROCESSORS, --processors PROCESSORS
                          threads (default: 6)
    -d, --debug           make extra debugging output (default: False)
    -h, --help            show this help message and exit

  PLOTTING:
    -pl [PLOTS [PLOTS ...]], --plots [PLOTS [PLOTS ...]]
                          Plots. Input 'all' or 'a' to plot all
                          1) Primary clustering dendrogram
                          2) Secondary clustering dendrograms
                          3) Secondary clustering MDS
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

  $ dRep evaluate -h
  usage: dRep evaluate [-p PROCESSORS] [-d] [-h] [--warn_dist WARN_DIST]
                     [--warn_sim WARN_SIM] [--warn_aln WARN_ALN]
                     [-e [EVALUATE [EVALUATE ...]]]
                     work_directory

  positional arguments:
  work_directory        Directory where data and output
                        *** USE THE SAME WORK DIRECTORY FOR ALL DREP OPERATIONS ***

  SYSTEM PARAMETERS:
  -p PROCESSORS, --processors PROCESSORS
                        threads (default: 6)
  -d, --debug           make extra debugging output (default: False)
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

Bonus
-----

Bonus consists of operations that don't really fit in with the functions of dRep, but can be helpful. Currently the only thing it can do is determine taxonomy of your bins. This is done using centrifuge, similar to how `anvi'o does it <http://merenlab.org/2016/06/18/importing-taxonomy/>`_. If you choose to use this option, the taxonomy of genome will be shown with the filename in most figures.
