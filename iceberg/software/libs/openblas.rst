.. Warning:: 
    Iceberg reaches end-of-life on **30th November 2020.**
    If you are running jobs on Iceberg then you need to take urgent action to ensure that your jobs/scripts will run on ShARC or Bessemer. 
 
    If you have never used ShARC or Bessemer then now is the time to test your scripts.
    Not all software on Iceberg is available on ShARC/Bessemer. 

.. _openblas_iceberg:

OpenBLAS
========

.. sidebar:: OpenBLAS

   :Latest version: 0.2.19
   :URL: http://www.openblas.net/

OpenBLAS is a library that provides low-level C and Fortran routines for linear algebra.  It provides optimised versions of the routines described in the `BLAS <https://en.wikipedia.org/wiki/Basic_Linear_Algebra_Subprograms>`_ specification.  It is based on GotoBLAS2 1.13 BSD version.

Usage
-----
To make OpenBLAS's shared libraries, static libraries and header files available you need to run: ::

        module load libs/gcc/4.9.2/openblas/0.2.19

**Note**: this will also load GCC 4.9.2 as OpenBLAS depends on the ``libgfortran``, ``libquadmath`` and ``libgcc_s`` libraries provided by that compiler.

You can then statically link with ``libopenblas.a`` or with ``-lopenblas`` if you want to use the shared library.

The library is **multithreaded**.  You can the number of threads with using one of the following environment variables: ::

* ``export OPENBLAS_NUM_THREADS=4``
* ``export GOTO_NUM_THREADS=4``
* ``export OMP_NUM_THREADS=4``

The priorities are ``OPENBLAS_NUM_THREADS`` > ``GOTO_NUM_THREADS`` > ``OMP_NUM_THREADS.``

**Note**: The build process console output indicates that the maximum number of threads is 12.

The library was compiled on and optimised for the Intel Nehalem architecture (inc. Intel Xeon X5650 CPUs) but should also run on the newer Intel Sandy Bridge-EP archtecture (Intel Xeon E5 2650/2660/2670 CPUs).

Documentation
-------------

See the project website including the `FAQ <https://github.com/xianyi/OpenBLAS/wiki/Faq>`_.

Installation notes
------------------
This section is primarily for administrators of the system. 

Version 0.2.19
^^^^^^^^^^^^^^

This was compiled with GCC 4.9.2 on an Intel Xeon X5650 CPU.

#. First, download, configure, build, test and install using :download:`this script </iceberg/software/install_scripts/libs/gcc/4.9.2/openblas/0.2.19/install.sh>`.
#. The default makefile target runs a test suite.  All tests should pass (see the :download:`build process console output </iceberg/software/install_scripts/libs/gcc/4.9.2/openblas/0.2.19/install.log>`.
#. Next, install :download:`this modulefile </iceberg/software/modulefiles/libs/gcc/4.9.2/openblas/0.2.19>` as ``/usr/local/modulefiles/libs/gcc/4.9.2/openblas/0.2.19`` 
