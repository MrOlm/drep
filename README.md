# drep
De-replication of microbial genomes

# Normal use case

drep.py filter workD -g genomelist

drep.py cluster workD

drep.py choose workD

drep.py analyze workD -cviz a

drep.py adjust workD something

# Module summaries

filter

# Dataframes

Bdb

:   Holds sequence locations and filenames
:   genome, location

Mdb

:   Holds raw information about MASH comparisons
:   genome1,genome2,dist,p,kmers,similarity

Ndb

:   Holds raw information about ANIn comparions
:   reference, querry, alignment_length, ani, reference_coverage, querry_coverage, alignment_coverage, MASH_cluster, ...