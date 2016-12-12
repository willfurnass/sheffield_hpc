.. _gromacs:

GROMACS
=======

.. sidebar:: GROMACS

   :Latest version: 2016.1
   :Dependancies: CUDA (optional)
   :URL: http://www.gromacs.org/


Interactive Usage
-----------------

After connecting to iceberg (see :ref:`ssh`),  start an interactive session with the ``qsh`` or ``qrsh`` command. To make GROMACS available in this session, run one of the following commands: ::

        module load apps/gcc/4.9.2/gromacs/2016.1-cuda-8.0.44
        module load apps/gcc/4.9.2/gromacs/2016.1-serial
        module load apps/gcc/4.9.2/gromacs/2016.1-serial-mkl

        module load apps/intel/15/gromacs/5.1.2
        module load apps/intel/15/gromacs/5.1.2-cuda-7.5.18

Installation notes
-------------------

Version 2016.1
^^^^^^^^^^^^^^

There are two different builds of this version, both of which were built using :ref:`GCC 4.9.2 <gcc_iceberg>` and :ref:`Boost 1.60.0 <boost>`.

 * **2016.1-serial**: Uses (internally built, single precision) FFTW3 for fast fourier transforms and uses OS-provided BLAS/LAPACK routines.  NB cannot build using FFTW3 for FFT and MKL for BLAS/LAPACK.  No GPU or MPI support.
 * **2016.1-cuda-8.0.44**: Uses CUDA 8.0.44.  May also use an (internally built, single precision) FFTW3 for fast fourier transforms and use an OS-provided BLAS/LAPACK routines.  No MPI support.

Both versions have GCC 4.9.2 as a run-time dependency.  The CUDA build also has CUDA 8 as a run-time dependency.

GROMACS can be built with OpenMPI support but this has not yet been done for this version.

Version 5.1
^^^^^^^^^^^

Compilation Choices:

* Use latest intel compilers
* ``-DGMX_MPI=on``
* ``-DCMAKE_PREFIX_PATH=/usr/local/packages6/apps/intel/15/gromcas``
* ``-DGMX_FFT_LIBRARY="fftw3"``
* ``-DGMX_BUILD_OWN_FFTW=ON``

The script used to build gromacs can be found :download:`here
</iceberg/software/install_scripts/apps/gromacs/install_gromacs_5.1.sh>`.
