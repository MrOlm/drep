Overview
========

dRep is a python program which performs rapid pair-wise comparison of genome sets (and other bells and whistles). It's primary purpose is for genome de-replication, but it can do a lot more.

Genome de-replication
---------------------

De-replication is the process of identifying groups of genomes that are the "same" in a genome set, and removing all but the "best" genome from each redundant set. How similar genomes need to be to be considered "same", how to determine which genome is "best", and other important decisions are discussed in :doc:`choosing_parameters`)

A common use for genome de-replication is the case of individual assembly of metagenomic data. If metagenomic samples are collected in a series (as is commonly done), a common way to assemble the short reads it with a "co-assembly". That is combining the reads from all samples and assembling them together. The problem with this is assembling similar strains together can severely fragment assemblies, precluding recovery of a good genome bin. An alternative option is to assemble each sample separately, and then "de-replicate" the bins from each assembly to make a final genome set.

.. image:: images/Figure0.jpg
