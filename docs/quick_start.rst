Quick Start
===========

The functionality of dRep is broken up into modules. The modules can be run separately (see :doc:`module_descriptions`), or together in workflows. To see a list of the available modules, check the help::

 $ dRep -h

                 ...::: dRep v2.4.0 :::...

  Matt Olm. MIT License. Banfield Lab, UC Berkeley. 2017

  Choose one of the operations below for more detailed help. See https://drep.readthedocs.io/en/latest/index.html for documentation
  Example: dRep dereplicate -h

  Workflows:
  compare         -> Perform rapid pair-wise comparison on a list of genomes
  dereplicate     -> De-replicate a list of genomes

  Single operations:
  filter          -> Filter a genome list based on size, completeness, and/or contamination
  cluster         -> Compare and cluster a genome list based on MASH and ANIn/gANI
  choose          -> Choose the best genome from each genome cluster
  evaluate        -> Evaluate genome de-replication
  bonus           -> Other random operations (determine taxonomy / check dependencies)

De-replication
---------------

De-replication is the process of identifying groups of genomes that are the "same" in a genome set, and removing all but the "best" genome from each redundant set. How similar genomes need to be to be considered "same", how the "best" genome is chosen,  and other options can be adjusted (see :doc:`choosing_parameters`)

To de-replicate a set of genomes, run the following command::

 $ dRep dereplicate outout_directory -g path/to/genomes/*.fasta

This will automatically de-replicate the genome list and produce lots of information about it.

.. seealso::
  :doc:`example_output`
    to view example output
  :doc:`choosing_parameters`
    for guidance changing parameters


Genome comparison
-----------------

dRep is able to perform rapid genome comparisons for a group of genomes and visualize their relatedness. For example::

 $ dRep compare output_directory -g path/to/genomes/*.fasta

For help understanding the output, see :doc:`example_output`

To change the comparison parameters, see :doc:`choosing_parameters`

.. seealso::
  :doc:`example_output`
    to view example output
  :doc:`choosing_parameters`
    for guidance changing parameters
