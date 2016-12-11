# dRep
De-replication of microbial genomes.

When time-series metagenomes are available,
often scientists do a co-assembly of all time-points together and then
bin those genomes using the time-series signal. Due to strain variation,
however, often a genome will assemble better in an individual sample than the
co-assembly. Often in the end what scientists want is a single set of genomes,
and so assembling each time-point seperately creates the problem of having to
"de-replicate" the genomes assembled at each time-point to the single, best
genome available from all assemblies. Drep can help the user identity cases
where the same genotype was assembled multiple times, and choose the best
genome for downstream analysis.

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

# Dependencies

Not all dependencies are required for all operations. Excecutable must be in the
system path.

MASH - https://genomebiology.biomedcentral.com/articles/10.1186/s13059-016-0997-x

Nucmer - http://mummer.sourceforge.net/

CheckM - http://ecogenomics.github.io/CheckM/

gANI - https://ani.jgi-psf.org/html/download.php?

centrifuge (bonus) - https://omictools.com/centrifuge-tool

# Normal use case

```
dRep dereplicate_wf workD -g genomelist
```

# DataTables

## The following datatables will be available in the work directory in the folder "data_tables" upon completion of a run

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

# Other user-facing data

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
