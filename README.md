# drep
De-replication of microbial genomes

# Normal use case

```
drep.py filter workD -g genomelist

drep.py cluster workD

drep.py choose workD

drep.py analyze workD -cviz a

drep.py adjust workD something
```

# Module summaries

## filter

Filters a genome list, Bdb, or Wdb, based on length and checkM information

Can modify:  
- Bdb, Chdb, Wdb

Can read:  
- Bdb, Wdb  

## cluster

Clusters genomes first by MASH, and then by a secondary method (ANIn or gANI)

Can modify:

- Cdb, Mdb, Ndb

Can read:

-  Bdb

## choose

Chooses the winning representative of each secondary cluster

Can modify:

- Bdb, Cdb, Chdb, Sdb, Wdb

Can read:

- Chdb 

## analyze

Make plots and visualize how new clustering methods would look

Can modify:

- NA

Can read:

-

# DataTables

Bdb

:   Holds sequence locations and filenames

:   genome, location

Mdb

:   Holds raw information about MASH comparisons

:   genome1,genome2,dist,p,kmers,similarity

Ndb

:   Holds raw information about ANIn comparions

:   reference, querry, alignment_length, ani, reference_coverage, querry_coverage, alignment_coverage, MASH_cluster, ...

# Other user-facing data

./figures

./dereplicated_genomes

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