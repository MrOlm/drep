Installation
============

To install dRep, simply run ::

$ pip install drep

OR ::

  $ git clone git@github.com:MrOlm/drep.git

  $ cd drep

  $ pip install .

That's it!

Dependencies
------------

dRep requires many other great programs to run. Not all dependencies are needed for all operations

* `Mash <https://genomebiology.biomedcentral.com/articles/10.1186/s13059-016-0997-x>`_
* `Nucmer <http://mummer.sourceforge.net/>`_
* `CheckM <http://ecogenomics.github.io/CheckM/>`_
* `gANI (aka CalculateANI) <https://ani.jgi-psf.org/html/download.php?>`_
* `Centrifuge <https://omictools.com/centrifuge-tool>`_

Testing
-------

To make sure everything is installed correctly you can run the dRep test suite using py.test::

 mattolm@biotite ~/Programs/drep $ cd drep/tests

 mattolm@biotite ~/Programs/drep/tests $ py.test

 === test session starts ===
 platform linux -- Python 3.5.1, pytest-3.0.5, py-1.4.32, pluggy-0.4.0
 rootdir: /home/mattolm/Programs/drep, inifile:
 collected 3 items

 test_suite.py ...

 === 3 passed in 506.46 seconds ===

pyenv
-----

Because dRep is written in python3 and CheckM is written in python2, you may need to use `pyenv <https://github.com/yyuu/pyenv>`_ to be able to call both.

With CheckM installed in a python2 installation of pyenv, and dRep installed in the python3 version, the following command should set allow both python2 and python3 commands to be called::

 $ pyenv global anaconda3-4.1.0 anaconda2-4.1.0
