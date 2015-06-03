iceberg Documentation
=====================

This repository builds the iceberg documentation using sphinx.

Building the .pdf file on Windows
---------------------------------
Install the following:-

* [Anaconda Python](https://store.continuum.io/cshop/anaconda/) 

Install the following module

     pip install sphinx_bootstrap_theme

* [GNU Make](http://gnuwin32.sourceforge.net/packages/make.htm)
* [MikTeX](http://miktex.org/download)

From the command line:

    make latexpdf

On first run, MikTeX will prompt you to install various extra LaTeX packages.

