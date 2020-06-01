.. _sharc-python-conda:

Python
======

.. sidebar:: Python

   :Support Level: Core
   :URL: https://python.org


This page documents the *Miniconda* installation on ShARC. This is the
recommended way of using Python, and the best way to be able to configure custom
sets of packages for your use.

*Conda* a Python package manager, allows you to create *environments* which are
sets of packages that you can modify. It does this by installing them in your
home area. This page will guide you through loading Conda and then creating and
modifying environments so you can install and use whatever Python packages you
need.

Using Conda Python
------------------

After connecting to ShARC (see :ref:`ssh`),
:ref:`start an interactive session <sched_interactive>` then
load Conda using: ::

   module load apps/python/conda

The ``root`` Conda environment (the default) provides Python 3 and no extra modules.
It is not recommended for general use, just as a base for your own environments.
There is also a ``python2`` environment,
which is the same but with a Python 2 installation.

Quickly Loading Anaconda Environments
-------------------------------------

There are a small number of environments provided for everyone to use.
These are the default ``root`` and ``python2`` environments
as well as various versions of Anaconda for Python 3 and Python 2.

The Anaconda environments can be loaded through provided module files: ::

   module load apps/python/anaconda2-4.2.0
   module load apps/python/anaconda3-4.2.0

Where ``anaconda2`` represents Python 2 installations and 
``anaconda3`` represents Python 3 installations.
These commands load the ``apps/python/conda`` module and then
activate the Anaconda environment specified.

.. note::
   Anaconda 2.5.0 and higher are compiled with Intel MKL libraries which should
   result in higher numerical performance.

Using Conda Environments
########################

Once the Conda module is loaded 
you have to load or create the desired Conda environments.
For the documentation on Conda environments see `the Conda documentation <http://conda.pydata.org/docs/using/envs.html>`_.

You can load a Conda environment with: ::

   source activate python2

where ``python2`` is the name of the environment, and unload one with::

   source deactivate

which will return you to the ``root`` environment.

It is possible to list all the available environments with: ::

   conda env list

Provided system-wide are a set of Anaconda environments.
These are installed with the Anaconda version number in the environment name, 
and are never modified.
They will therefore provide a static base for derivative environments or for using directly.

Creating an Environment
#######################

Every user can create their own environments, and packages shared with the
system-wide environments will not be reinstalled or copied to your file store: 
they will be ``symlinked``, which reduces the space you need in your ``/home``
directory to install many different Python environments.

To create a clean environment with just Python 2 and numpy you can run: ::

   conda create -n mynumpy python=2.7 numpy

This will download the latest release of Python 2.7 and numpy, and create an
environment named ``mynumpy``.

Any version of Python or list of packages can be provided: ::

   conda create -n myscience python=3.5 numpy=1.8.1 scipy

If you wish to modify an existing environment, such as one of the Anaconda
installations, you can ``clone`` that environment: ::

   conda create --clone anaconda3-4.2.0 -n myexperiment

This will create an environment called ``myexperiment`` which has all the
Anaconda 4.2.0 packages installed with Python 3.


Installing Packages Inside an Environment
#########################################

Once you have created your own environment 
you can install additional packages or different versions of packages into it.
There are two methods for doing this: ``conda`` and ``pip``
If a package is available through Conda it is
strongly recommended that you use Conda to install that package.
You can search for packages using Conda: ::

   conda search pandas

then install the package using: ::

   conda install pandas

if you are not in your environment
you will get a permission denied error when trying to install packages.
If this happens, create or activate an environment you own.

If a package is not available through Conda you can search for and install it
using pip e.g.: ::

   pip search colormath

   pip install colormath


Using Python with MPI
---------------------

There is an **experimental** set of packages for Conda
that have been compiled by the RSE and RCG teams,
which allow you to use a MPI stack entirely managed by Conda.
This allows you to easily create complex environments and 
use MPI without worrying about other modules or system libraries.

To get access to these packages you need to 
run the following command to add the repo to your Conda config: ::

   conda config --add channels file:///usr/local/packages/apps/conda/conda-bld/

you should then be able to install the packages with the ``openmpi`` feature,
which currently include ``openmpi``, ``hdf5``, ``mpi4py`` and ``h5py``: ::

   conda create -n my_mpi_env python=3.5 openmpi mpi4py

Currently, this channel provides Conda packages for:

 - ``mpi4py`` (and ``openmpi``) for Python 3.4, 3.5, 3.6 and 2.7
 - ``h5py`` (and ``hdf5``) with MPI support for Python 3.5 and 2.7

The build scripts for these packages can be found in 
this `GitHub <https://github.com/rcgsheffield/conda-packages>`_ repository.

Installation Notes
------------------

These are primarily for administrators of the system.

The Conda package manager is installed in ``/usr/share/packages/apps/conda``. 
It was installed using the `miniconda <http://conda.pydata.org/miniconda.html>`_ installer.

It is important to regularly update the ``root`` environment 
to keep the Conda package manager up to date. 
To do, ensure you have write permissions on ``/usr/local/packages/apps/conda`` then run: ::

   $ conda update --all
   $ conda update conda

Installing a New Version of Anaconda
####################################

Run the following as a ``sa_`` user (with write permissions to
``/usr/local/packages/apps/conda``: ::

   $ conda create -n anaconda3-<VERSION> python=3 anaconda=<VERSION>
   $ conda create -n anaconda2-<VERSION> python=2 anaconda=<VERSION>

Then copy the modulefile for the previous version of Anaconda to the new
version and update the name of the environment. Also you will need to append
the new module to the ``conflict`` line in `apps/python/.conda-environments.tcl`.
