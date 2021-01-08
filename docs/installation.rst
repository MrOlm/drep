Installation
============

Using pip
---------

To install dRep, simply run ::

$ pip install drep

OR ::

  $ git clone https://github.com/MrOlm/drep.git

  $ cd drep

  $ pip install .

That's it!

Pip is a great package with many options to change the installation parameters in various ways. For details, see `pip documentation <https://packaging.python.org/installing/>`_

Using conda
----------------

To install dRep with conda, simply run ::

  conda config --add channels bioconda; conda install drep

Dependencies
------------

dRep requires other programs to run. Not all dependencies are needed for all operations

To check which dependencies are installed on your system and accessible by dRep, run ::

 $ dRep check_dependencies

**Near Essential**

* `Mash <https://genomebiology.biomedcentral.com/articles/10.1186/s13059-016-0997-x>`_ - Makes primary clusters (v1.1.1 confirmed works)
* `MUMmer <http://mummer.sourceforge.net/>`_ - Performs default ANIm comparison method (v3.23 confirmed works)

**Recommended**

* `fastANI <https://github.com/ParBLiSS/FastANI>`_ - A fast secondary clustering algorithm
* `CheckM <http://ecogenomics.github.io/CheckM/>`_ - Determines contamination and completeness of genomes (v1.0.7 confirmed works)
* `gANI (aka ANIcalculator) <https://ani.jgi-psf.org/html/download.php?>`_ - Performs gANI comparison method (v1.0 confirmed works)
* `Prodigal <http://prodigal.ornl.gov/>`_ - Used be both checkM and gANI (v2.6.3 confirmed works)

**Accessory**

* `NSimScan <https://pubmed.ncbi.nlm.nih.gov/27153714/>`_ - Only needed for goANI algorithm
* `Centrifuge <https://omictools.com/centrifuge-tool>`_ - Deprecated; not used anymore

Programs need to be installed to the system path, so that you can call them from anywhere on your computer.

.. note::

  If you already have information on your genome's completeness and contamination, you can input that to dRep without the need to install checkM (see :doc:`advanced_use`))


CheckM
-------

CheckM is the program that dRep uses to calculate completeness and contamination. It's not required for all dRep commands, but if you'd like to de-dereplicate your genome set using completeness and contamination it is required. It takes a bit of work to install (including setting a root directory and downloading a reference database). Installation instructions can be found `here <https://github.com/Ecogenomics/CheckM/wiki/Installation#how-to-install-checkm>`_
