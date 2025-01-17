.. _fftw_stanage:

fftw
====

.. sidebar:: fftw

   :Latest version: 3.3.10
   :URL: https://www.fftw.org/

FFTW is a C subroutine library for
computing the discrete Fourier transform (DFT)
in one or more dimensions,
of arbitrary input size,
and of both real and complex data
(as well as of even/odd data,
i.e. the discrete cosine/sine transforms or DCT/DST).

Usage
-----

.. include:: /referenceinfo/imports/scheduler/SLURM/common_commands/srun_start_interactive_session_import_stanage.rst

To make this library available, run one the following:

.. include:: /referenceinfo/imports/stanage/packages/fftw-ml-el7-icelake-znver-stanage.rst

- `gompi` versions are a subset of the :ref:`foss toolchain <stanage_eb_toolchains>`
  and also load GCC and OpenMPI. :ref:`See matching foss toolchains<foss-toolchain-table>`.
- `gompic` versions are a subset of the :ref:`fosscuda toolchain <stanage_eb_toolchains>`
  and also load GCC, OpenMPI and CUDA.


Also see :ref:`imkl-fftw <imkl_fftw_stanage>` which is a library that combines FFTW library with Intel's Math Kernel Library (IMKL)
to provide optimized FFT routines that are specifically optimized for Intel processors.

.. include:: /referenceinfo/imports/stanage/packages/fftw-dpnd-el7-icelake-znver-stanage.rst
