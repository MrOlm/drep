# dRep

[![Downloads](https://pepy.tech/badge/drep)](https://pepy.tech/project/drep)
[![Downloads](https://pepy.tech/badge/drep/week)](https://pepy.tech/project/drep)

dRep is a python program for rapidly comparing large numbers of genomes. dRep can also "de-replicate" a genome set by identifying groups of highly similar genomes and choosing the best representative genome for each genome set.

Manual, installation instructions, and API are at available at
[ReadTheDocs](https://drep.readthedocs.io/en/latest/)

Publication is available at
[ISMEJ](https://doi.org/10.1038/ismej.2017.126)

Open source pre-print publication is available at
[bioRxiv](https://doi.org/10.1101/108142)

## ⚡ New in v4

dRep v4 uses [skani](https://github.com/bluenote-1577/skani) for **both** primary and secondary genome comparisons by default, replacing v3's default of MASH (primary) + fastANI (secondary). skani is much faster than that pair, and it *streams* its comparisons instead of building an all-vs-all matrix in memory — so `dereplicate` runs far quicker and its memory footprint grows roughly linearly with genome count rather than quadratically (the O(N²) scaling that previously limited large runs — see [#259](https://github.com/MrOlm/drep/issues/259)).

**Whole-pipeline `dRep dereplicate` on 10,000 genomes** — identical inputs and settings:

| | v3 defaults (MASH → fastANI) | v4 defaults (skani) |
|---|---|---|
| Wall-clock time | 6 h 15 min | **14.5 min** (~26× faster) |
| Peak memory (RSS) | ~13 GB | ~7 GB |

*Benchmarked on an Apple M1 Pro with `-p 10`; the memory advantage widens further at larger genome counts.*

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
* [skani](https://github.com/bluenote-1577/skani) - Makes primary clusters and performs the default secondary comparison (v0.2+ confirmed works)
* [CheckM](http://ecogenomics.github.io/CheckM/) - Determines contamination and completeness of genomes (v1.0.7 confirmed works). Only needed for `dereplicate`; skip it with `--genomeInfo` or `--ignoreGenomeQuality`

### Optional

* [Mash](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-016-0997-x>) - Only needed for `--primary_algorithm MASH` (v1.1.1 confirmed works)
* [MUMmer](http://mummer.sourceforge.net/) - Only needed for the ANIm comparison methods (v3.23 confirmed works)
* [fastANI](https://github.com/ParBLiSS/FastANI) - An alternative fast secondary clustering algorithm
* [gANI (aka ANIcalculator)](https://ani.jgi-psf.org/html/download.php?) - Performs gANI comparison method (v1.0 confirmed works)
* [Prodigal](http://prodigal.ornl.gov/) - Used be both checkM and gANI (v2.6.3 confirmed works)
* [NSimScan](https://pubmed.ncbi.nlm.nih.gov/27153714/) - Only needed for goANI algorithm (open source version of gANI)
