Quick Start
===========

The functionality of dRep is broken up into modules. To see a list of the available modules, check the help::

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


Dereplication
---------------

Dereplication is the process of identifying groups of genomes that are the "same" in a genome set and identifying the "best" genome from each set. How similar genomes need to be to be considered "same" and how the "best" genome is chosen are study-specific decisions that can be adjusted (see :doc:`choosing_parameters`)

To dereplicate a set of genomes, run the following command::

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
