.. image:: https://travis-ci.org/rcgsheffield/sheffield_hpc.svg?branch=master
    :target: https://travis-ci.org/rcgsheffield/sheffield_hpc
.. image:: https://readthedocs.org/projects/iceberg/badge/?version=latest
    :target: https://readthedocs.org/projects/iceberg/builds/

Sheffield High Performance Computing Documentation
==================================================

This is the source code for the documentation of Sharc and Iceberg, The University of Sheffield's High Performance Computing clusters. It is written in the rst format.

For a guide on the rst file format see `this <http://thomas-cokelaer.info/tutorials/sphinx/rest_syntax.html>`_ document.

Rendered Documentation
----------------------
Two versions of the documentation are currently automatically built from this repository:

* `A website <http://docs.hpc.shef.ac.uk/en/latest/>`_.
* `A .pdf document <http://readthedocs.org/projects/iceberg/downloads/pdf/latest/>`_.

How to Contribute
-----------------
To contribute to this documentation, first you have to fork it on GitHub and clone it to your machine, see `Fork a Repo <https://help.github.com/articles/fork-a-repo/>`_ for the GitHub documentation on this process.

Once you have the git repository locally on your computer, you will need to install ``sphinx`` and ``sphinx_bootstrap_theme`` to be able to build the documentation. See the instructions below for how to achieve this.

Once you have made your changes and updated your Fork on GitHub you will need to `Open a Pull Request <https://help.github.com/articles/using-pull-requests/>`_.
All changes to the repository should be made through Pull Requests, including those made by the people with direct push access.


Building the documentation on a local Windows machine
#####################################################

#. Install **Python 3.6** on your machine by downloading and running the `Miniconda for Python 3.6 <https://conda.io/miniconda.html>`_ installer: 

    * Install for *just you*;
    * Install to the default location (e.g. ``C:\Users\myusername\Miniconda3``);
    * Do **not** *add Anaconda to your PATH environment variable*;
    * Do **not** *register Anaconda as your default Python 3.6*.

#. Click *Start* -> *Anaconda3 (64-bit)* -> *Anaconda Prompt* to open a terminal window.

#. Create a new *conda environment* for building the documentation by running the following from this window: ::

    conda create --name sheffield_hpc python=3.6
    conda activate sheffield_hpc	# . activate sheffield_hpc on older versions of conda
    pip install -r requirements.txt

#. To build the HTML documentation run: ::

    make html
	
   Or if you don't have the ``make`` utility installed on your machine then build with *sphinx* directly: ::

    sphinx-build . ./html

#. If you want to build the PDF documentation you will need:

    * `GNU Make <http://gnuwin32.sourceforge.net/packages/make.htm>`_
    * `MikTeX <http://miktex.org/download>`_

   Then from the command line, the following will build the ``.pdf`` file: ::

    make latexpdf

   On first run, MikTeX will prompt you to install various extra LaTeX packages.

Building the documentation on a local Linux machine
###################################################

#. Ensure Python 3 (ideally Python 3.6) is installed.
#. Create a `virtual environment <https://docs.python.org/3/tutorial/venv.html>`_ to install sphinx into: ::

    mkdir -m 700 ~/.venvs
    python3 -m venv ~/.venvs/sheffield_hpc_py3
    source ~/.venvs/sheffield_hpc_py3/bin/activate

#. Install the Python packages needed to build the HTML documentation: ::

     pip3 install -r requirements.txt

#. Build the documentation: ::

     make html

Building the documentation on a local Mac machine
#################################################

#. Ensure Python 3 (ideally Python 3.6) is installed.  If you do not already have a python distribution installed, we recommend you install `Miniconda for Python 3.6 <https://conda.io/miniconda.html>`_.
#. Install the Python packages needed to build the HTML documentation.  If you are using (mini)conda create a new *conda environment* for building the documentation by running: ::

    export PATH=${HOME}/miniconda3/bin:$PATH
    conda create -n sheffield_hpc python=3.6
    pip install -r requirements.txt

   If you are *not* using (mini)conda to provide Python 3: ::

    mkdir -m 700 ~/.venvs
    python3 -m venv ~/.venvs/sheffield_hpc_py3
    source ~/.venvs/sheffield_hpc_py3/bin/activate
    pip3 install --requirement requirements.txt

#. To build the HTML documentation run::

    make html

Continuous build and serve
##########################

The package `sphinx-autobuild <https://github.com/GaretJax/sphinx-autobuild>`_ provides a watcher that automatically rebuilds the site as files are modified. To use it, install (in addition to the Sphinx packages) with the following: ::

    pip install sphinx-autobuild

To start the autobuild process, run: ::

    make livehtml

The application also serves up the site at port ``8000`` by default at http://localhost:8000.


Making Changes to the Documentation
-----------------------------------

The documentation consists of a series of `reStructured Text <http://sphinx-doc.org/rest.html>`_ files which have the ``.rst`` extension.
These files are then automatically converted to HTML and combined into the web version of the documentation by sphinx.
It is important that when editing the files the syntax of the rst files is followed.
If there are any errors in your changes the build will fail and the documentaion  will not update, you can test your build locally by running ``make html``.
The easiest way to learn what files should look like is to read the ``rst`` files already in the repository.
