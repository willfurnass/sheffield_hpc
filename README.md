Iceberg Documentation
=====================

Iceberg is The University of Sheffield's High Performance Computing cluster. The current official documentation for it is at [https://www.shef.ac.uk/wrgrid/iceberg](https://www.shef.ac.uk/wrgrid/iceberg).

This repository contains a proposed replacement for the page linked to above. It makes use of [Sphinx](http://sphinx-doc.org/) to build the results and represents a move to a more open, collaborative way of working.

Rendered Documentation
----------------------
Two versions of the documentation are currently automatically built from this repository:

* (A website)[http://rcg.group.shef.ac.uk/iceberg/]
* (A .pdf document)[http://rcg.group.shef.ac.uk/iceberg/icebergDocumentation.pdf]

How to Contribute
-----------------
Contribution is via githib Pull Requests, even for those who have direct commit access to the repository. To contribute, please fork this repo, make your changes in the fork and submit a Pull Request when you are ready.

Building the documentation locally
----------------------------------

It is also possible to build the documentation on your own machine. Clone the repository and follow the instructions below

Building the .pdf file on a local Windows machine
-------------------------------------------------

Install the following:-

* [Anaconda Python](https://store.continuum.io/cshop/anaconda/) 

Install the following module

     pip install sphinx_bootstrap_theme

* [GNU Make](http://gnuwin32.sourceforge.net/packages/make.htm)
* [MikTeX](http://miktex.org/download)

From the command line:

    make latexpdf

On first run, MikTeX will prompt you to install various extra LaTeX packages.

