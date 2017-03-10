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

Dependencies
------------

dRep requires many other great programs to run. Not all dependencies are needed for all operations

* `Mash <https://genomebiology.biomedcentral.com/articles/10.1186/s13059-016-0997-x>`_
* `Nucmer <http://mummer.sourceforge.net/>`_
* `CheckM <http://ecogenomics.github.io/CheckM/>`_
* `Prodigal <http://prodigal.ornl.gov/>`_
* `gANI (aka CalculateANI) <https://ani.jgi-psf.org/html/download.php?>`_
* `Centrifuge <https://omictools.com/centrifuge-tool>`_

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
