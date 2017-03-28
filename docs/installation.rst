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

Dependencies
------------

dRep requires many other great programs to run. Not all dependencies are needed for all operations

**Near Essential**

* `Mash <https://genomebiology.biomedcentral.com/articles/10.1186/s13059-016-0997-x>`_ - Makes primary clusters (v1.1.1 confirmed works)
* `MUMmer <http://mummer.sourceforge.net/>`_ - Performs ANIm comparison method (v3.23 confirmed works)

**Recommended**

* `CheckM <http://ecogenomics.github.io/CheckM/>`_ - Determines contamination and completeness of genomes (v1.0.7 confirmed works)
* `gANI (aka CalculateANI) <https://ani.jgi-psf.org/html/download.php?>`_ - Performs gANI comparison method (v1.0 confirmed works)
* `Prodigal <http://prodigal.ornl.gov/>`_ - Used be both checkM and gANI (v2.6.3 confirmed works)

**Accessory**

* `Centrifuge <https://omictools.com/centrifuge-tool>`_ - Performs taxonomic assignment of bins (v1.0.3 confirmed works)

Programs need to be installed to the system path, so that you can call them from anywhere on your computer.

Testing
-------

To make sure everything is installed correctly you can run the dRep test suite::

 $ cd drep/tests

 $ python test_suite.py

pyenv
-----

Because dRep is written in python3 and CheckM is written in python2, you may need to use `pyenv <https://github.com/yyuu/pyenv>`_ to be able to call both.

With CheckM installed in a python2 installation of pyenv, and dRep installed in the python3 version, the following command should set allow both python2 and python3 commands to be called::

 $ pyenv global 3.5.1 2.7.9

Alternatively, you could add python2 to your CheckM shebang line (though I have not confirmed that this works)
