.. _theano_iceberg:

Theano
======

.. sidebar:: Theano

   :URL: http://deeplearning.net/software/theano/index.html

Theano is a Python library that allows you to define, optimize, and evaluate mathematical expressions involving multi-dimensional arrays efficiently.
Theano is most commonly used to perform deep learning and has excellent GPU support and integration through PyCUDA.
The following steps can be used to setup and configure Theano on your own profile.

About Theano on Iceberg
-----------------------

As Theano and all its dependencies are written in Python, it can be installed locally in your home directory.
The use of Anaconda (:ref:`python-conda`) is recommended as it is able to create a virtual environment in your home directory,
allowing the installation of new Python packages without admin permission.

Installation
------------

First :ref:`start an interactive session <sched_interactive>`.
To use GPUs see :ref:`GPUInteractive_iceberg`.

Load the relevant modules with the following commands: ::

   module load apps/python/anaconda3-2.5.0
   module load libs/binlibs/cudnn/5.1-cuda-8.0.44

Create a Conda environment to load relevant modules on your local user account: ::

   conda create -n theano python=3.5
   source activate theano

Install the other Python module dependencies which are required using pip (alternatively these could be installed with Conda if you prefer) ::

   pip install theano nose nose-parameterized pycuda

For optimal Theano performance, enable the CUDA memory manager *CNMeM*.
To do this, create/update an ``.theanorc`` file in your home directory and
set the fraction of GPU memory reserved by Theano.
The exact amount of energy may have to be hand-picked:
if Theano asks for more memory that is currently available on the GPU,
an error will be thrown during import of theano module.

For example to reserve 80% of GPU memory your ``.theanorc`` would contain: ::

   [lib]
   cnmem=0.8

Run python and verify that Theano is working correctly: ::

   python -c "import theano;theano.test()"

Every Session Afterwards and in Your Job Scripts
------------------------------------------------

The instructions above install Theano and its dependencies inside your home directory but
every time you use a new session or within your job scripts,
the modules must be loaded and Conda must be activated again.
Use the following command to activate the Conda environment with Theano installed: ::

   module load apps/python/anaconda3-2.5.0
   module load libs/binlibs/cudnn/5.1-cuda-8.0.44
   source activate theano
