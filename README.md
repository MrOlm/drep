# dRep
De-replication of microbial genomes.

Genome-resolved metagenomic studies often target groups of samples from the same environment to resolve ecosystem processes and to provide constraints for series-based binning. Many studies co-assemble reads from all samples to avoid redundancies that arise because the same genome is reconstructed from multiple samples. Here we present dRep, a program that identifies redundancies and selects the highest quality genome from the redundant set. By sequentially applying a fast inaccurate estimation of genome distance and a slow but accurate measure of average nucleotide identity, the program reduces the computational time for de-replication by several orders of magnitude. The assembly of data from individual samples followed by de-replication increases the number of high-quality genomes reconstructed and decreases their fragmentation relative to co-assembly. dRep is an open-source python program freely available under a BSD license.

Pre-print publication available at bioRxiv:
[link to publication](https://www.youtube.com/watch?v=dQw4w9WgXcQ)

# Installation

```
pip install drep
```

OR

```
git clone git@github.com:MrOlm/drep.git

cd drep

pip install .
```

## Dependencies

**Not all dependencies are required for all operations.** Executable must be in the
system path.

- [Mash](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-016-0997-x)

- [Nucmer](http://mummer.sourceforge.net/)

- [CheckM](http://ecogenomics.github.io/CheckM/)

- [gANI](https://ani.jgi-psf.org/html/download.php?)

- [centrifuge](https://omictools.com/centrifuge-tool)

## Testing

To make sure everything is installed correctly you can run the dRep test suite using py.test

```
mattolm@biotite ~/Programs/drep $ cd drep/tests

mattolm@biotite ~/Programs/drep/tests $ py.test

=== test session starts ===
platform linux -- Python 3.5.1, pytest-3.0.5, py-1.4.32, pluggy-0.4.0
rootdir: /home/mattolm/Programs/drep, inifile:
collected 3 items

test_suite.py ...

=== 3 passed in 506.46 seconds ===

```

## pyenv

Because dRep is written in python3 and CheckM is written in python2, you may need to use [pyenv](https://github.com/yyuu/pyenv) to be able to call both.

With CheckM installed in a python2 installation of pyenv, and dRep installed in the python3 version, the following command should set allow both python2 and python3 commands to be called:

```
pyenv global anaconda3-4.1.0 anaconda2-4.1.0
```

# Quick start

There are several steps that dRep uses in order to de-replicate a genome set. The user can perform each step individually (see Advanced use), or run a "workflow" to run multiple steps together with one command. For a list of possible commands, run the help:

```
mattolm@biotite /home $ dRep -h

                ...::: dRep v0.3.0 :::...

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
```

## Genome De-replication

To de-replicate a list of genomes using default settings, use dereplicate_wf:

```
dRep dereplicate_wf test_output -g path/to/genomes/*.fasta
```

To see the list of parameters, run:

```
dRep dereplicate_wf -h
```

## Comparing genomes without de-replication

To simply compare a list of genomes using Mash and a secondary algorithm, use compare_wf:

```
dRep compare_wf test_output -g path/to/genomes/*.fasta
```

# Basic use

## Steps of de-replication

There are several steps that dRep uses in order to de-replicate a genome set. The user can perform each step individually, or run a "workflow" to run multiple steps together with one command (see Quick Start).

- **filter**: Filter the genome list based on size, completeness, and/or contamination

- **cluster**: Compare and cluster the genome list based on both Mash and a secondary algorithm (ANIm / gANI)

- **choose**: Choose the best genome from each cluster based on CheckM

- **analyze**: Make figures based on the de-replication process above

There are  all additional steps that can be performed in addition to those listed above, including determining taxonomy of bins (**bonus**), evaluating genome de-replication (**evaluate**), and adjusting the thresholds for specific clusters (**adjust**).

## The work directory

This is a little bit of information about the work directory

## Interpreting output

This is the basics of what the output is and how to interpret it

### DataTables

The following datatables will be available in the work directory in the folder "data_tables" upon completion of a run

Bdb

:   Holds sequence locations and filenames

:   genome, location

Mdb

:   Holds raw information about MASH comparisons

:   genome1, genome2, dist, p, kmers, similarity

Ndb

:   Holds raw information about ANIn comparions

:   reference, querry, alignment_length, ani, reference_coverage, querry_coverage, alignment_coverage, MASH_cluster, ...

Cdb

:   Holds clustering information

:   cluster_method, comparison_algorithm, genome, primary_cluster, secondary_cluster, threshold

### Other user-facing data

./figures

./dereplicated_genomes

./Clustering_files

These pickle files store information on both primary and secondary clusters.
Loading the **first** value gives you the linkage, loading the **second** value gives you the db that was used to make the linkage, loading the **third** value give you a dictionary of the arguments that were used to make the linkage.

Example:  
```
f = open(pickle, 'rb')
linkage = pickle.load(f)
db = pickle.load(f)
arguments = pickle.load(f)
```

# Work directory tree

workDirectory  
./data  
...../MASH_files/  
...../ANIn_files/  
...../gANI_files/  
...../Clustering_files/  
...../checkM/  
........./genomes/  
........./checkM_outdir/  
...../prodigal/  
./figures  
./data_tables  
...../Bdb.csv  # Sequence locations and filenames  
...../Mdb.csv  # Raw results of MASH comparisons  
...../Ndb.csv  # Raw results of ANIn comparisons  
...../Cdb.csv  # Genomes and cluster designations  
...../Chdb.csv # CheckM results for Bdb  
...../Sdb.csv  # Scoring information  
...../Wdb.csv  # Winning genomes  
./dereplicated_genomes  
./log  
...../logger.log  
...../cluster_arguments.json  

# Known bugs

* gANI cannot handle comparisons where scaffolds have the same names. Be sure that
all scaffold names in all genomes are unique BEFORE the first non-fasta character
(space, colon, etc.)
