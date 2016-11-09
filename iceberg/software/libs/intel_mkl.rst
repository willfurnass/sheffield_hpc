.. _iceberg_intel_mkl:

Intel Math Kernel Library
=========================

Intel's Math Kernel Library (MKL) provides highly optimized, threaded and
vectorized functions to maximize performance on each processor family. It
Utilises de-facto standard C and Fortran APIs for compatibility with BLAS,
LAPACK and FFTW functions from other math libraries.

Parallel Studio Composer Edition version
----------------------------------------

The MKL was installed along with the :ref:`Intel compilers
<iceberg_intel_compilers>` as part of Intel Parallel Studio.  The MKL can be
used with or without other Parallel Studio packages (such as the Intel
compilers).  

To activate just the MKL, use one of : ::

    module load libs/binlibs/intel-mkl/2017.0

    module load libs/binlibs/intel-mkl/11.2.3

Note that 2017.0 is newer than 11.2.3.

Installation Notes
------------------

The following notes are primarily for system administrators.

**Intel MKL 2017.0**

Installed as part of :ref:`Parallel Studio Composer Edition 2017
<iceberg_intel_parallel_studio>`.

:download:`This modulefile
</iceberg/software/modulefiles/libs/binlibs/intel-mkl/2017.0>` was installed as
``/usr/local/modulefiles/libs/binlibs/intel-mkl/2017.0``.

**Intel MKL 11.2.3**

Installed as part of :ref:`Intel Parallel Studio Composer Edition 2015 Update 3
<iceberg_intel_compilers>`.

:download:`This modulefile
</iceberg/software/modulefiles/libs/binlibs/intel-mkl/11.2.3>` was installed as
``/usr/local/modulefiles/libs/binlibs/intel-mkl/11.2.3``
