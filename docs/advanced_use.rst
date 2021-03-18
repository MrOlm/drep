Advanced Use
============

Accessing Internal Information
------------------------------

All of internal information is stored in the work directory

work directory file-tree
+++++++++++++++++++++++++

::

  workDirectory
  ./data
  ...../checkM/
  ...../Clustering_files/
  ...../gANI_files/
  ...../MASH_files/
  ...../ANIn_files/
  ...../prodigal/
  ./data_tables
  ...../Bdb.csv  # Sequence locations and filenames
  ...../Cdb.csv  # Genomes and cluster designations
  ...../Chdb.csv # CheckM results for Bdb
  ...../Mdb.csv  # Raw results of MASH comparisons
  ...../Ndb.csv  # Raw results of ANIn comparisons
  ...../Sdb.csv  # Scoring information
  ...../Wdb.csv  # Winning genomes
  ...../Widb.csv # Winning genomes' checkM information
  ./dereplicated_genomes
  ./figures
  ./log
  ...../cluster_arguments.json
  ...../logger.log
  ...../warnings.txt

Data Tables
+++++++++++

Within the ``data_tables`` folder is where organized data lives. It's very easy to access this information, as it's all stored in .csv files.

.. note::
  If you code in Python, I cannot recommend `pandas <http://pandas.pydata.org/>`_ enough for data-frame reading and manipulation. It's how all data is manipulated behind the scenes in dRep. See the API section below for easy access to these dataframes

Bdb
  Genome input locations, filenames, and lengths
Cdb
  Primary cluster, Secondary cluster, and information on clustering method for each genome
Chdb
  CheckM results for all genomes
Mdb
  Pair-wise Mash comparison results
Ndb
  Secondary comparison results
Tdb
  Taxonomy (as determined by centrifuge)
Sdb
  The score of each genome
Wdb
  The cluster and score of de-replicated genomes
Widb
  Useful checkM information on de-replicated genomes

Clustering files
++++++++++++++++

These pickle files store information on both primary and secondary clusters. Loading the first value gives you the linkage, loading the second value gives you the db that was used to make the linkage, loading the third value give you a dictionary of the arguments that were used to make the linkage.

For example::

  f = open(pickle, 'rb')
  linkage = pickle.load(f)
  db = pickle.load(f)
  arguments = pickle.load(f)

Raw data
++++++++

Refer to the above file structure to find the rest of the raw data. The data is kept from all program runs, so you can find the raw ANIm/gANI files, Mash sketches, prodigal gene predictions, centrifuge raw output, ect.

Using external genome quality information
--------

If you already have your own genome quality information and would not like dRep to run checkM to generate it again, you can provide it using the `genomeInformation` flag.

The genomeInformation file must be in .csv format and have the columns "genome", "completeness", and "contamination". Columns "completeness" and "contamination" should be 0-100, and "genome" is the filename of the genome.

For example::

  genome,completeness,contamination
  Enterococcus_casseliflavus_EC20.fasta,98.28,0.0
  Enterococcus_faecalis_T2.fna,98.28,0.0
  Enterococcus_faecalis_TX0104.fa,96.55,0.0
  Enterococcus_faecalis_YI6-1.fna,98.28,0.0
  Escherichia_coli_Sakai.fna,100.0,0.0

Caching
--------

The reason that dRep stores all of the raw information that it does is so that if future operations require the same file, it can just load the file that's already there instead of making it again. This goes for prodigal gene predictions, checkM, centrifuge, all ANI algorithms, ect. The data-frame that uses the information **will** be remade, but the information itself will not.

The reason I mention this is because if you would like to run another dRep operation that's similar to one you've already run, you can use the results of the last run to speed up the second run.

For example, say you've already run the dereplicate_wf using gANI and want to run the same thing but with ANIm to compare. You can make a copy of the gANI work directory, and then run the dereplicate_wf on the copy specifying the new secondary algorithm. It will have to run all of the ANIm comparisons, but will **not** re-run checkM, prodigal, centrifuge, ect., as the data will already be cached in the work directory.

.. warning::

  Be warned, this is somewhat buggy and can easily get out of hand. While it does save time, sometimes it's just best to re-run the whole thing again with a clean start

Restarting a clustering job with already completed primary clustering
++++++

There are (rare) circumstances where you may way to run primary clustering on one machine and secondary clustering on another. Primary cluster needs more RAM, and secondary clustering needs more cores, so doing this could let you choose the optimal machine for each specific step. There's not a formal way of making dRep do this, but there is an hacky way. If you run dRep in debug mode (with ``-d``), it'll make a file named "CdbF.csv" in the data_tables folder. This file contains the primary clustering information, and has a very simple format that looks like this::

    genome,primary_cluster
    Enterococcus_casseliflavus_EC20.fasta,0
    Enterococcus_faecalis_T2.fna,0
    Enterococcus_faecalis_TX0104.fa,0
    Enterococcus_faecalis_YI6-1.fna,0
    Escherichia_coli_Sakai.fna,1

If you run a dRep job on a work directort that already has this file in the ``data_tables`` folder, as well as a file named ``Mdb.csv`` (it doesn't actually matter what's in that file), as well as in debug mode (``-d``), dRep will load this file instead of running primary clustering. You can also make the CdbF.csv file yourself using the format above.

I know this is confusing- it's a hacky thing I threw together for my own research, but I figured I would let everyone else know about it as well. Feel free to shoot me an email if this is something you're interested in doing and can't get it to work.

API
---

See :doc:`source/drep` for the API to dRep. For example::

  from drep.WorkDirectory import WorkDirectory

  wd = WorkDirectory('path/to/workdirectory')
  Mdb = wd.get_db('Mdb')
  Cdb = wd.get_db('Cdb')
  ...

This will work for all datatables.

Be warned that the API is not very well maintained or documented, and you'll likely have to do a bit of digging into the source code if you want to use it extensively.

Troubleshooting checkM
------------------------

One of the most common problems users have when running dRep are failures related to the program checkM. These errors can show up as either of the following::

  New checkM db needs to be made

  !!! checkM failed !!!

These errors can be caused by checkM crashing a variety of ways; here are some tips to figure out what's wrong.

.. note::

  Sometimes the easier thing to do is just run checkM (or whatever genome assessment tool you prefer) yourself and provide the results to dRep, instead of making dRep run checkM. See the above section "Using external genome quality information" for more info.

1) Ensure that you have checkM installed. You can do this with the `check_dependencies option`::

    With dRep version 3+

    $ dRep check_dependencies
    mash.................................... all good        (location = /home/mattolm/miniconda3/envs/drep_testing/bin/mash)
    nucmer.................................. all good        (location = /home/mattolm/miniconda3/envs/drep_testing/bin/nucmer)
    checkm.................................. all good        (location = /home/mattolm/miniconda3/envs/drep_testing/bin/checkm)
    ANIcalculator........................... !!! ERROR !!!   (location = None)
    prodigal................................ all good        (location = /home/mattolm/miniconda3/envs/drep_testing/bin/prodigal)
    centrifuge.............................. all good        (location = /home/mattolm/miniconda3/envs/drep_testing/bin/centrifuge)
    nsimscan................................ !!! ERROR !!!   (location = None)
    fastANI................................. all good        (location = /home/mattolm/bin/fastANI)`

    With dRep version 2+

    $ dRep bonus test --check_dependencies
    mash.................................... all good        (location = /home/mattolm/miniconda3/envs/drep_testing/bin/mash)
    nucmer.................................. all good        (location = /home/mattolm/miniconda3/envs/drep_testing/bin/nucmer)
    checkm.................................. all good        (location = /home/mattolm/miniconda3/envs/drep_testing/bin/checkm)
    ANIcalculator........................... !!! ERROR !!!   (location = None)
    prodigal................................ all good        (location = /home/mattolm/miniconda3/envs/drep_testing/bin/prodigal)
    centrifuge.............................. all good        (location = /home/mattolm/miniconda3/envs/drep_testing/bin/centrifuge)
    nsimscan................................ !!! ERROR !!!   (location = None)
    fastANI................................. all good        (location = /home/mattolm/bin/fastANI)`

If checkM reports ``!! ERROR !!``, try and re-install it.

2) Make sure you have have the checkM data installed and accessible, as described here: https://github.com/Ecogenomics/CheckM/wiki/Installation#how-to-install-checkm

3) To see the specific error that checkM is throwing, re-run dRep with the ``-d`` flag. This will produce a number of files in the ``log/cmd_logs/`` folder. Looking through these will be the commands that dRep gave to checkM, and the STDERR and STDOUT that checkM produced. Looking at the actual error code checkM is giving can be really helpful in figuring out what's wrong

4) If you're running lots of genomes, sometimes python will hit a recursion limit while running checkM and stall out. To fix this you can increase the recursion limit by setting the flag ``--set_recursion`` to some really big number. **Note: as of dRep version 3.2.0, dRep runs checkM in groups to prevent this problem. See ``checkm_group_size`` for more info.**

5) A newer problem is the error ``New checkM db needs to be made``. This usually means that checkM worked on some, but not all of your genomes. Check ``log/logger.log`` to see which ones failed. I haven't quite figured out this problem yet. If you encounter it I would encourage you to just run checkM outside of dRep. If you think you know how to fix it, please send me an email.

6) You can check checkM's run log at ''data/checkM/checkM_outdir/''. Sometimes there's helpful info there as to what's going wrong.

7) Check dRep's log file for lines with the words ``Running CheckM with command``. This will tell you the exact checkM command dRep tried to run. Trying to run these commands yourself can be another way of seeing what the actual error is and troubleshooting.

8) If all else fails, please post a GitHub issue. I'm happy to help troubleshoot.


