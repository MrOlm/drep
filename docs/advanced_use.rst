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

Caching
--------

The reason that dRep stores all of the raw information that it does is so that if future operations require the same file, it can just load the file that's already there instead of making it again. This goes for prodigal gene predictions, checkM, centrifuge, all ANI algorithms, ect. The data-frame that uses the information **will** be remade, but the information itself will not.

The reason I mention this is because if you would like to run another dRep operation that's similar to one you've already run, you can use the results of the last run to speed up the second run.

For example, say you've already run the dereplicate_wf using gANI and want to run the same thing but with ANIm to compare. You can make a copy of the gANI work directory, and then run the dereplicate_wf on the copy specifying the new secondary algorithm. It will have to run all of the ANIm comparisons, but will **not** re-run checkM, prodigal, centrifuge, ect., as the data will already be cached in the work directory.

.. warning::

  Be warned, this is somewhat buggy and can easily get out of hand. While it does save time, sometimes it's just best to re-run the whole thing again with a clean start

API
---

This is no formal API for dRep, but of course you can read through the source-code and import any function that you'd like. The one class that is particularly helpful, however, is the WorkDirectory. By loading the results of your dRep run as a WorkDirectory object, you have direct access to all of the data-frames.

For example::

  from drep.WorkDirectory import WorkDirectory

  wd = WorkDirectory('path/to/workdirectory')
  Mdb = wd.get_db('Mdb')
  Cdb = wd.get_db('Cdb')
  ...

This will work for all datatables
