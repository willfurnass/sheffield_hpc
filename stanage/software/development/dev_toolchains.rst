.. _stanage_eb_toolchains:

Development toolchains
======================

Much of the software available on Stanage
was installed using `EasyBuild <https://easybuild.readthedocs.io/>`__,
which is a tool for building and installing software on HPC systems.

EasyBuild groups core sets of dependencies as what it calls *toolchains*.
The two supported toolchains on Stanage are:

- ``foss``
   - C, C++ and Fortran compilers: :ref:`GCC <gcc_stanage>`
   - MPI implementation: :ref:`OpenMPI <openmpi_stanage>`
   - BLAS and LAPACK implementation: :ref:`OpenBLAS <openblas_stanage>`
   - Parallel, distributed LAPACK implementation: :ref:`ScaLAPACK <scalapack_stanage>`
   - Fourier transforms: :ref:`FFTW <fftw_stanage>`

- ``intel``
   - C, C++ and Fortran compilers (:ref:`icc/icpc/ifort <iccifort_stanage>`)
   - MPI implementation (:ref:`Intel MPI <impi_stanage>`)
   - BLAS, LAPACK and fourier transforms: :ref:`Intel MKL <imkl_stanage>`

Sub-toolchains
--------------

In addition to the above toolchains there are sub-toolchains 
corresponding to subsets of the main toolchain dependencies:

* ``gompi``: ``GCC`` + ``openmpi``
* ``gompic``: ``GCC`` + ``OpenMPI`` + ``CUDA``
* ``gcccuda``: ``GCC`` + ``CUDA``
* ``iccifort``: ``icc`` + ``ifort``
* ``iimpi``: ``icc`` + ``ifort`` + ``impi``

Dependency versions for toolchains
----------------------------------

See the `EasyBuild documentation for dependency versions for foss and intel <https://docs.easybuild.io/common-toolchains>`__

.. _foss-toolchain-table:

.. table:: 
   :widths: 25 25 25 25

   ============== ====== =========  ====== 
   Foss toolchain GCC    ScaLAPACK  FFTW   
   ============== ====== =========  ====== 
   2024a          13.3.0 2.2.0      3.3.10   
   2023b          13.2.0 2.2.0      3.3.10 
   2023a          12.3.0 2.2.0      3.3.10 
   2022b          12.2.0 2.2.0      3.3.10 
   2022a          11.3.0 2.2.0      3.3.10 
   2021b          11.2.0 2.1.0      3.3.10 
   2021a          10.3.0 2.1.0      3.3.9
   ============== ====== =========  ====== 
   
.. _intel-toolchain-table:

.. table:: 
   :widths: 25 25 25 25

   =============== ====== ========= =========
   Intel toolchain  GCC   IMPI      IMKL   
   =============== ====== ========= =========
   2024a           13.3.0 2021.13.0 2024.2.0 
   2023b           13.2.0 2021.10.1 2023.2.0       
   2023a           12.3.0 2021.9.1  2023.1.0         
   2022b           12.2.0 2021.7.1  2022.2.1                
   2022a           11.3.0 2021.6.1  2022.1.0                
   2021b           11.2.0 2021.4.0  2021.4.0                
   2021a           10.3.0 2021.2.0  2021.2.0                
   =============== ====== ========= =========

