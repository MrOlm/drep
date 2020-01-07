# dRep

dRep is a python program for rapidly comparing large numbers of genomes. dRep can also "de-replicate" a genome set by identifying groups of highly similar genomes and choosing the best representative genome for each genome set.

Manual, installation instructions, and API are at available at
[ReadTheDocs](https://drep.readthedocs.io/en/latest/)

Publication is available at
[ISMEJ](http://www.nature.com/ismej/journal/vaop/ncurrent/full/ismej2017126a.html)

Open source pre-print publication is available at
[bioRxiv](https://doi.org/10.1101/108142)

## Installation with pip
```
$ pip install drep
```

## Quick start

### Genome comparison:
```
$ dRep compare output_directory -g path/to/genomes/*.fasta
```

### Genome de-replication:
```
$ dRep dereplicate outout_directory -g path/to/genomes/*.fasta
```

### Make sure dependencies are properly installed:
```
$ drep bonus output_directory --check_dependencies
```

## Dependencies
### Required
* [Mash](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-016-0997-x) is used to rapidly compare all genomes in a pair-wise manner
* [MUMmer](http://mummer.sourceforge.net/) is used to perform more actuate comparisons between genomes which are shown to be similar with Mash

### Optional
* [CheckM](http://ecogenomics.github.io/CheckM/) is used to determine the contamination and completeness of genomes (used during de-replication)
* [gANI (aka ANIcalculator)](https://ani.jgi-psf.org/html/download.php?) is an optional alternative to MUMmer
* [Prodigal](http://prodigal.ornl.gov/) is a dependency of both checkM and gANI

### Accessory
* [Centrifuge](https://omictools.com/centrifuge-tool) can be used to perform rough taxonomic assignment of bins
