User manual
===================

dRep has 3 commands: compare, dereplicate, and check dependencies. To see a list of these options check the help::

  $ dRep -h

                ...::: dRep v3.0.0 :::...

  Matt Olm. MIT License. Banfield Lab, UC Berkeley. 2017 (last updated 2020)

  See https://drep.readthedocs.io/en/latest/index.html for documentation
  Choose one of the operations below for more detailed help.

  Example: dRep dereplicate -h

  Commands:
    compare            -> Compare and cluster a set of genomes
    dereplicate        -> De-replicate a set of genomes
    check_dependencies -> Check which dependencies are properly installed

In previous versions of dRep (everything before v3) the user could run a number of additional modules separately, but now they can only be run as part of the larger workflows `compare` and `dereplicate`. Many of the modules are the same for `compare` and `dereplicate`, however, and in cases where these is the same parameter in both it functions exactly the same in each.

dRep has descriptions in the program help for all the adjustable parameters. If any of these are particularly confusing, don't hesitate to send an email to ask what it does.


.. seealso::

  :doc:`Important Concepts`
    for theoretical thoughts about how to choose appropriate parameters and thresholds

  :doc:`example_output`
     for help interpreting the output from your run in the work directory

  :doc:`advanced_use`
    for access to the raw output data and the python API

Compare
--------

This workflow compares a set of genomes. For a list of all parameters, check the help::

  $ dRep compare -h
    usage: dRep compare [-p PROCESSORS] [-d] [-h] [-g [GENOMES [GENOMES ...]]]
                        [--S_algorithm {fastANI,gANI,goANI,ANIn,ANImf}]
                        [-ms MASH_SKETCH] [--SkipMash] [--SkipSecondary]
                        [--n_PRESET {normal,tight}] [-pa P_ANI] [-sa S_ANI]
                        [-nc COV_THRESH] [-cm {total,larger}]
                        [--clusterAlg {median,weighted,single,complete,average,ward,centroid}]
                        [--multiround_primary_clustering]
                        [--primary_chunksize PRIMARY_CHUNKSIZE]
                        [--greedy_secondary_clustering]
                        [--run_tertiary_clustering] [--warn_dist WARN_DIST]
                        [--warn_sim WARN_SIM] [--warn_aln WARN_ALN]
                        work_directory

    positional arguments:
      work_directory        Directory where data and output are stored
                            *** USE THE SAME WORK DIRECTORY FOR ALL DREP OPERATIONS ***

    SYSTEM PARAMETERS:
      -p PROCESSORS, --processors PROCESSORS
                            threads (default: 6)
      -d, --debug           make extra debugging output (default: False)
      -h, --help            show this help message and exit

    GENOME INPUT:
      -g [GENOMES [GENOMES ...]], --genomes [GENOMES [GENOMES ...]]
                            genomes to filter in .fasta format. Not necessary if
                            Bdb or Wdb already exist. Can also input a text file
                            with paths to genomes, which results in fewer OS
                            issues than wildcard expansion (default: None)

    GENOME COMPARISON OPTIONS:
      --S_algorithm {fastANI,gANI,goANI,ANIn,ANImf}
                            Algorithm for secondary clustering comaprisons:
                            fastANI = Kmer-based approach; very fast
                            ANImf   = (DEFAULT) Align whole genomes with nucmer; filter alignment; compare aligned regions
                            ANIn    = Align whole genomes with nucmer; compare aligned regions
                            gANI    = Identify and align ORFs; compare aligned ORFS
                            goANI   = Open source version of gANI; requires nsmimscan
                             (default: ANImf)
      -ms MASH_SKETCH, --MASH_sketch MASH_SKETCH
                            MASH sketch size (default: 1000)
      --SkipMash            Skip MASH clustering, just do secondary clustering on
                            all genomes (default: False)
      --SkipSecondary       Skip secondary clustering, just perform MASH
                            clustering (default: False)
      --n_PRESET {normal,tight}
                            Presets to pass to nucmer
                            tight   = only align highly conserved regions
                            normal  = default ANIn parameters (default: normal)

    GENOME CLUSTERING OPTIONS:
      -pa P_ANI, --P_ani P_ANI
                            ANI threshold to form primary (MASH) clusters
                            (default: 0.9)
      -sa S_ANI, --S_ani S_ANI
                            ANI threshold to form secondary clusters (default:
                            0.95)
      -nc COV_THRESH, --cov_thresh COV_THRESH
                            Minmum level of overlap between genomes when doing
                            secondary comparisons (default: 0.1)
      -cm {total,larger}, --coverage_method {total,larger}
                            Method to calculate coverage of an alignment
                            (for ANIn/ANImf only; gANI and fastANI can only do larger method)
                            total   = 2*(aligned length) / (sum of total genome lengths)
                            larger  = max((aligned length / genome 1), (aligned_length / genome2))
                             (default: larger)
      --clusterAlg {median,weighted,single,complete,average,ward,centroid}
                            Algorithm used to cluster genomes (passed to
                            scipy.cluster.hierarchy.linkage (default: average)

    GREEDY CLUSTERING OPTIONS
    These decrease RAM use and runtime at the expense of a minor loss in accuracy.
    Recommended when clustering 5000+ genomes:
      --multiround_primary_clustering
                            Cluster each primary clunk separately and merge at the
                            end with single linkage. Decreases RAM usage and
                            increases speed, and the cost of a minor loss in
                            precision and the inability to plot
                            primary_clustering_dendrograms. Especially helpful
                            when clustering 5000+ genomes. Will be done with
                            single linkage clustering (default: False)
      --primary_chunksize PRIMARY_CHUNKSIZE
                            Impacts multiround_primary_clustering. If you have
                            more than this many genomes, process them in chunks of
                            this size. (default: 5000)
      --low_ram_primary_clustering
                            Use a memory-efficient algorithm for primary clustering.
                            This only affects primary clustering and not secondary
                            clustering. Can be combined with multiround_primary_clustering
                            for even greater memory efficiency. (default: False)
      --greedy_secondary_clustering
                            Use a heuristic to avoid pair-wise comparisons when
                            doing secondary clustering. Will be done with single
                            linkage clustering. Only works for fastANI S_algorithm
                            option at the moment (default: False)
      --run_tertiary_clustering
                            Run an additional round of clustering on the final
                            genome set. This is especially useful when greedy
                            clustering is performed and/or to handle cases where
                            similar genomes end up in different primary clusters.
                            Only works with dereplicate, not compare. (default:
                            False)

    WARNINGS:
      --warn_dist WARN_DIST
                            How far from the threshold to throw cluster warnings
                            (default: 0.25)
      --warn_sim WARN_SIM   Similarity threshold for warnings between dereplicated
                            genomes (default: 0.98)
      --warn_aln WARN_ALN   Minimum aligned fraction for warnings between
                            dereplicated genomes (ANIn) (default: 0.25)

    Example: dRep compare output_dir/ -g /path/to/genomes/*.fasta


Dereplicate
------------

This workflow dereplicates a set of genomes. For a list of all parameters, check the help::

    $ dRep dereplicate -h
    usage: dRep dereplicate [-p PROCESSORS] [-d] [-h] [-g [GENOMES [GENOMES ...]]]
                            [-l LENGTH] [-comp COMPLETENESS] [-con CONTAMINATION]
                            [--ignoreGenomeQuality] [--genomeInfo GENOMEINFO]
                            [--checkM_method {taxonomy_wf,lineage_wf}]
                            [--set_recursion SET_RECURSION]
                            [--S_algorithm {goANI,ANIn,gANI,ANImf,fastANI}]
                            [-ms MASH_SKETCH] [--SkipMash] [--SkipSecondary]
                            [--n_PRESET {normal,tight}] [-pa P_ANI] [-sa S_ANI]
                            [-nc COV_THRESH] [-cm {total,larger}]
                            [--clusterAlg {single,ward,complete,weighted,centroid,median,average}]
                            [--multiround_primary_clustering]
                            [--primary_chunksize PRIMARY_CHUNKSIZE]
                            [--greedy_secondary_clustering]
                            [--run_tertiary_clustering]
                            [-comW COMPLETENESS_WEIGHT]
                            [-conW CONTAMINATION_WEIGHT]
                            [-strW STRAIN_HETEROGENEITY_WEIGHT] [-N50W N50_WEIGHT]
                            [-sizeW SIZE_WEIGHT] [-centW CENTRALITY_WEIGHT]
                            [--warn_dist WARN_DIST] [--warn_sim WARN_SIM]
                            [--warn_aln WARN_ALN]
                            work_directory

    positional arguments:
      work_directory        Directory where data and output are stored
                            *** USE THE SAME WORK DIRECTORY FOR ALL DREP OPERATIONS ***

    SYSTEM PARAMETERS:
      -p PROCESSORS, --processors PROCESSORS
                            threads (default: 6)
      -d, --debug           make extra debugging output (default: False)
      -h, --help            show this help message and exit

    GENOME INPUT:
      -g [GENOMES [GENOMES ...]], --genomes [GENOMES [GENOMES ...]]
                            genomes to filter in .fasta format. Not necessary if
                            Bdb or Wdb already exist. Can also input a text file
                            with paths to genomes, which results in fewer OS
                            issues than wildcard expansion (default: None)

    GENOME FILTERING OPTIONS:
      -l LENGTH, --length LENGTH
                            Minimum genome length (default: 50000)
      -comp COMPLETENESS, --completeness COMPLETENESS
                            Minumum genome completeness (default: 75)
      -con CONTAMINATION, --contamination CONTAMINATION
                            Maximum genome contamination (default: 25)

    GENOME QUALITY ASSESSMENT OPTIONS:
      --ignoreGenomeQuality
                            Don't run checkM or do any quality filtering. NOT
                            RECOMMENDED! This is useful for use with
                            bacteriophages or eukaryotes or things where checkM
                            scoring does not work. Will only choose genomes based
                            on length and N50 (default: False)
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

    GENOME COMPARISON OPTIONS:
      --S_algorithm {goANI,ANIn,gANI,ANImf,fastANI}
                            Algorithm for secondary clustering comaprisons:
                            fastANI = Kmer-based approach; very fast
                            ANImf   = (DEFAULT) Align whole genomes with nucmer; filter alignment; compare aligned regions
                            ANIn    = Align whole genomes with nucmer; compare aligned regions
                            gANI    = Identify and align ORFs; compare aligned ORFS
                            goANI   = Open source version of gANI; requires nsmimscan
                             (default: ANImf)
      -ms MASH_SKETCH, --MASH_sketch MASH_SKETCH
                            MASH sketch size (default: 1000)
      --SkipMash            Skip MASH clustering, just do secondary clustering on
                            all genomes (default: False)
      --SkipSecondary       Skip secondary clustering, just perform MASH
                            clustering (default: False)
      --n_PRESET {normal,tight}
                            Presets to pass to nucmer
                            tight   = only align highly conserved regions
                            normal  = default ANIn parameters (default: normal)

    GENOME CLUSTERING OPTIONS:
      -pa P_ANI, --P_ani P_ANI
                            ANI threshold to form primary (MASH) clusters
                            (default: 0.9)
      -sa S_ANI, --S_ani S_ANI
                            ANI threshold to form secondary clusters (default:
                            0.95)
      -nc COV_THRESH, --cov_thresh COV_THRESH
                            Minmum level of overlap between genomes when doing
                            secondary comparisons (default: 0.1)
      -cm {total,larger}, --coverage_method {total,larger}
                            Method to calculate coverage of an alignment
                            (for ANIn/ANImf only; gANI and fastANI can only do larger method)
                            total   = 2*(aligned length) / (sum of total genome lengths)
                            larger  = max((aligned length / genome 1), (aligned_length / genome2))
                             (default: larger)
      --clusterAlg {single,ward,complete,weighted,centroid,median,average}
                            Algorithm used to cluster genomes (passed to
                            scipy.cluster.hierarchy.linkage (default: average)

    GREEDY CLUSTERING OPTIONS
    These decrease RAM use and runtime at the expense of a minor loss in accuracy.
    Recommended when clustering 5000+ genomes:
      --multiround_primary_clustering
                            Cluster each primary clunk separately and merge at the
                            end with single linkage. Decreases RAM usage and
                            increases speed, and the cost of a minor loss in
                            precision and the inability to plot
                            primary_clustering_dendrograms. Especially helpful
                            when clustering 5000+ genomes. Will be done with
                            single linkage clustering (default: False)
      --primary_chunksize PRIMARY_CHUNKSIZE
                            Impacts multiround_primary_clustering. If you have
                            more than this many genomes, process them in chunks of
                            this size. (default: 5000)
      --low_ram_primary_clustering
                            Use a memory-efficient algorithm for primary clustering.
                            This only affects primary clustering and not secondary
                            clustering. Can be combined with multiround_primary_clustering
                            for even greater memory efficiency. (default: False)
      --greedy_secondary_clustering
                            Use a heuristic to avoid pair-wise comparisons when
                            doing secondary clustering. Will be done with single
                            linkage clustering. Only works for fastANI S_algorithm
                            option at the moment (default: False)
      --run_tertiary_clustering
                            Run an additional round of clustering on the final
                            genome set. This is especially useful when greedy
                            clustering is performed and/or to handle cases where
                            similar genomes end up in different primary clusters.
                            Only works with dereplicate, not compare. (default:
                            False)

    SCORING CRITERIA
    Based off of the formula:
    A*Completeness - B*Contamination + C*(Contamination * (strain_heterogeneity/100)) + D*log(N50) + E*log(size) + F*(centrality - S_ani)

    A = completeness_weight; B = contamination_weight; C = strain_heterogeneity_weight; D = N50_weight; E = size_weight; F = cent_weight:
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
      -centW CENTRALITY_WEIGHT, --centrality_weight CENTRALITY_WEIGHT
                            Weight of (centrality - S_ani) (default: 1)

    WARNINGS:
      --warn_dist WARN_DIST
                            How far from the threshold to throw cluster warnings
                            (default: 0.25)
      --warn_sim WARN_SIM   Similarity threshold for warnings between dereplicated
                            genomes (default: 0.98)
      --warn_aln WARN_ALN   Minimum aligned fraction for warnings between
                            dereplicated genomes (ANIn) (default: 0.25)

    Example: dRep dereplicate output_dir/ -g /path/to/genomes/*.fasta

Work Directory
--------------

The work directory is where all of the program's internal workings, log files, cached data, and output is stored.

.. seealso::

  :doc:`example_output`
    for help finding where the output from your run is located in the work directory

  :doc:`advanced_use`
    for access to the raw internal data (which can be very useful)

Genome filtering
-----------------

In the `dereplicate` module, the genome set is quality filtered first (for why this is necessary, see :doc:`choosing_parameters`). This is done using checkM. All genomes which don't pass the length threshold are filtered first to avoid running checkM unnecessarily. All genomes which don't pass checkM thresholds are filtered before comparisons are run to avoid running comparisons unnecessarily.

.. warning::

  All genomes must have at least one ORF called or else checkM will stall, so a length minimum of at least 10,000bp is recommended.

Warnings
--------

A series of checks are preformed to alert the user to potential problems with de-replication. There are two things that it looks for:

**de-replicated genome similarity**- this is comparing all of the de-replicated genomes to each other and making sure they're not too similar. This is to try and catch cases where similar genomes were split into different primary clusters, and thus failed to be de-replicated. *Depending on the number of de-replicated genomes, this can take a while*

**secondary clusters that were almost different**- this alerts you to cases where genomes are on the edge between being considered "same" or "different", depending on the clustering parameters you used. *This module reads the parameters you used during clustering from the work directory, so you don't need to specify them again.*

Overall these warnings are a bit half-baked, however, and I personally don't pay attention to them when running dRep myself.