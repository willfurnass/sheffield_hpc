.. Warning:: 
    Iceberg reaches end-of-life on **30th November 2020.**
    If you are running jobs on Iceberg then you need to take urgent action to ensure that your jobs/scripts will run on ShARC or Bessemer. 
 
    If you have never used ShARC or Bessemer then now is the time to test your scripts.
    Not all software on Iceberg is available on ShARC/Bessemer. 

.. _`Intel R`:

R (Intel Build)
===============

.. sidebar:: R

   :Dependencies: BLAS
   :URL: http://www.r-project.org/
   :Documentation: http://www.r-project.org/

R is a statistical computing language. This version of R is built using the :ref:`Intel compilers <iceberg_intel_compilers>` and the Intel Math Kernel Library. This combination can result in significantly faster runtimes in certain circumstances.

Most R extensions are written and tested for the gcc suite of compilers so it is recommended that you perform testing before switching to this version of R.

CPU Architecture
----------------
The Intel build of R makes use of CPU instructions that are only present on the most modern of Iceberg's nodes. In order to use them, you should add the following to your batch submission scripts:

.. code-block:: bash

   #$ -l arch=intel-e5-2650v2

If you do not do this, you will receive the following error message when you try to run R: ::

    Please verify that both the operating system and the processor support Intel(R) F16C and AVX instructions.

Loading the modules
-------------------
After connecting to iceberg (see :ref:`ssh`),  start an interactive session with the :code:`qrshx` command.

There are two types of the Intel builds of R, ``sequential`` and ``parallel`.` ``sequential`` makes use of one CPU core and can be used as a drop-in replacement for the standard version of R installed on Iceberg.

.. code-block:: bash

   module load apps/intel/15/R/3.3.1_sequential

The ``parallel`` version makes use of multiple CPU cores for certain linear algebra routines since it is linked to the parallel version of the Intel MKL. Note that **only** linear algebra routines are automatically parallelised.

.. code-block:: bash

   module load apps/intel/15/R/3.3.1_parallel

When using the parallel module, you must also ensure that you set the bash environment variable ``OMP_NUM_THREADS`` to the number of cores you require and also use the openmp parallel environment.  E.g. Add the following to your submission script:

.. code-block:: bash

   #$ -pe openmp 8
   export OMP_NUM_THREADS=8

   module load apps/intel/15/R/3.3.1_parallel


Example batch jobs
------------------
Here is how to run the R script called ``linear_algebra_bench.r`` from the `HPC Examples <https://github.com/mikecroucher/HPC_Examples>`_ github repository:

.. code-block:: bash

   #!/bin/bash
   #This script runs the linear algebra benchmark multiple times using the intel-compiled version of R
   #that's linked with the sequential MKL
   #$ -l rmem=8G
   # Target the Ivy Bridge Processors
   #$ -l arch=intel-e5-2650v2
 
   module load apps/intel/15/R/3.3.1_sequential
 
   echo "Intel R with sequential MKL on intel-e5-2650v2"
   Rscript linear_algebra_bench.r output_data.rds

Here is how to run the same code using 8 cores:

.. code-block:: bash

   #!/bin/bash
   #$ -l rmem=3G # Memory per core
   # Target the Ivy Bridge Processors
   #$ -l arch=intel-e5-2650v2
   #$ -pe openmp 8
   export OMP_NUM_THREADS=8
 
   module load apps/intel/15/R/3.3.1_parallel
 
   echo "Intel R with parallel MKL on intel-e5-2650v2"
   echo "8 cores"
   Rscript inear_algebra_bench.r 8core_output_data.rds

Installing additional packages
------------------------------
By default, the standard version of R allows you to install packages into the location ``~/R/x86_64-unknown-linux-gnu-library/3.3/``, where ``~`` refers to your home directory.

To ensure that the Intel builds do not contaminate the standard gcc builds, the Intel R module files set the environment variable ``R_LIBS_USER`` to point to ``~/R/intel_R/3.3.1``

As a user, you should not need to worry about this detail and just install packages as you usuall would from within R. e.g.

.. code-block:: R

   install.packages("dplyr")

The Intel build of R will ignore any packages installed in your home directory for the standard version of R and vice versa

Installation Notes
------------------
These notes are primarily for administrators of the system.

version 3.3.1
^^^^^^^^^^^^^

This was a scripted install. It was compiled from source with Intel Compiler 15.0.3 and with ``--enable-R-shlib`` enabled. It was run in batch mode.

This build required several external modules including :ref:`xzutils`, :ref:`curl`, :ref:`bzip2` and :ref:`zlib`

Sequential
""""""""""

#. Install using :download:`this script </iceberg/software/install_scripts/apps/intel/15/R/sequential-3.3.1/install.sh>`. Downloads, compiles, tests and installs R 3.3.1 using Intel Compilers and the sequential MKL. The install and test logs are at ``/usr/local/packages6/apps/intel/15/R/sequential-3.3.1/install_logs/``
#. :download:`This module file </iceberg/software/modulefiles/apps/intel/15/R/3.3.1_sequential>` installed as ``apps/intel/15/R/3.3.1_sequential``

Parallel
""""""""

#. Install using :download:`this script </iceberg/software/install_scripts/apps/intel/15/R/parallel-3.3.1/install.sh>`.  Downloads, compiles, tests and installs R 3.3.1 using Intel Compilers and the parallel MKL. The install and test logs are at ``/usr/local/packages6/apps/intel/15/R/parallel-3.3.1/install_logs/`` 
#. :download:`This module file </iceberg/software/modulefiles/apps/intel/15/R/3.3.1_parallel>` installed as ``apps/intel/15/R/3.3.1_parallel``
