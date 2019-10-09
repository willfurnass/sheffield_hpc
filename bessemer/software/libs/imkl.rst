.. _imkl_bessemer:

Intel MKL
=========

.. sidebar:: Intel MKL
   
   :URL: https://software.intel.com/en-us/mkl
   :Documentation: https://software.intel.com/en-us/mkl/documentation/view-all

Intel's Math Kernel Library (MKL) provides
highly optimized, threaded and vectorized functions to
maximize performance on each processor family.
It utilises de-facto standard C and Fortran APIs
for compatibility with :ref:`BLAS <blas_bessemer>`,
LAPACK and
FFTW functions from other math libraries.

Usage
-----

The Intel MKL can be activated using: ::

   module load imkl/2019.1.144-iimpi-2019a

which also loads the :ref:`iimpi-2019a <bessemer_eb_toolchains>` toolchain.
