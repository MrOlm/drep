# dRep

[![Downloads](https://pepy.tech/badge/drep)](https://pepy.tech/project/drep)
[![Downloads](https://pepy.tech/badge/drep/week)](https://pepy.tech/project/drep)

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
$ dRep dereplicate output_directory -g path/to/genomes/*.fasta
```

### Make sure dependencies are properly installed:
```
$ dRep check_dependencies
```

## Dependencies
### Near Essential
* [Mash](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-016-0997-x>) - Makes primary clusters (v1.1.1 confirmed works)
* [MUMmer](http://mummer.sourceforge.net/) - Performs default ANIm comparison method (v3.23 confirmed works)

### Optional

* [fastANI](https://github.com/ParBLiSS/FastANI) - A fast secondary clustering algorithm
* [CheckM](http://ecogenomics.github.io/CheckM/)_ - Determines contamination and completeness of genomes (v1.0.7 confirmed works)
* [gANI (aka ANIcalculator)](https://ani.jgi-psf.org/html/download.php?) - Performs gANI comparison method (v1.0 confirmed works)
* [Prodigal](http://prodigal.ornl.gov/) - Used be both checkM and gANI (v2.6.3 confirmed works)
* [NSimScan](https://pubmed.ncbi.nlm.nih.gov/27153714/) - Only needed for goANI algorithm (open source version of gANI)
