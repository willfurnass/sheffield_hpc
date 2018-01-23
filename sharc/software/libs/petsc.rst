.. _petsc_sharc:

PETSc
=====

.. sidebar:: PETSc

   :Latest version: 3.8.3
   :URL: https://www.mcs.anl.gov/petsc/

PETSc, pronounced PET-see (the S is silent), is a suite of data structures and routines for the scalable (parallel) solution of scientific applications modeled by partial differential equations. 
It supports MPI, and GPUs through CUDA or OpenCL, as well as hybrid MPI-GPU parallelism. 
PETSc (sometimes called PETSc/Tao) also contains the Tao optimization software library. 

Note: the build installed on SHARC supports MPI but does not include CUDA/OpenCL support.

To get more of a general overview of PETSc see the PETSc website or listen to `this RCECast episode <http://www.rce-cast.com/Podcast/rce-24-petsc.html>`_.

Usage
-----

To make this library available, run one of the following: ::

   module load libs/petsc/3.8.3/intel-17.0.0-openmpi-2.0.1
   module load libs/petsc/3.8.3/gcc-6.2-openmpi-2.0.1

This loads a build of PETSc that has MPI support (but no CUDA support). 
It also loads 

* OpenMPI 2.0.1 built using the Intel compiler suite 17.0.0 **or**
* OpenMPI 2.0.1 built using the GCC compiler suite 6.2.

This build uses the Intel MKL that came with Intel Parallel Studio 2017.0 for BLAS/LAPACK functionality.

Examples and tutorials
----------------------

PETSc examples are saved in a particular location on SHARC: ::

   $PETSC_DIR/share/petsc/examples

and you want to copy them to a writeable location (e.g. :ref:`your /home or /data directory or $TMPDIR <fastdata>`) before compiling and running them. 
For example: ::

   cp -r $PETSC_DIR/share/petsc/examples /data/$USER/

Tutorials that use these examples are available `here <http://www.mcs.anl.gov/petsc/documentation/tutorials/HandsOnExercise.html>`_.
Note that these tutorials assume your examples are in a particular directory and don't tell you to copy the example files to a writable directory first. 

Here's a :ref:`batch job submission script <sge-queue>` for running one of the tutorial examples on ShARC using 16 CPU cores and MPI: ::

   #!/bin/bash
   #$ -pe mpi 16
   #$ -m bea 
   #$ -P insigneo-imsb
   #$ -q insigneo-imsb.q
   #$ -M some.user@sheffield.ac.uk
   #$ -j y

   # Load a particular version of PETSc (plus OpenMPI)
   module load libs/petsc/3.8.3/intel-17.0.0-openmpi-2.0.1

   # Copy the tutorials to a writable location
   cp -r $PETSC_DIR/share/petsc/examples /data/$USER/

   cd /data/$USER/examples/src/ts/examples/tutorials
   # Compile the example
   make ex2

   mpiexec -n $NSLOTS ./ex2 -ts_max_steps 10 -ts_monitor -M 128 

The `expected output <http://www.mcs.anl.gov/petsc/petsc-current/src/ts/examples/tutorials/output/ex2_tut_3.out.html>`_ for this example.

Note that not all examples may run as PETSc was not built with all possible solvers.

Installation notes
------------------
This section is primarily for administrators of the system. 

3.8.3 (Intel build)
^^^^^^^^^^^^^^^^^^^

* :download:`Install script</sharc/software/install_scripts/libs/petsc/3.8.3/intel-17.0.0-openmpi-2.0.1/install.sh>`;
* :download:`Install log</sharc/software/install_scripts/libs/petsc/3.8.3/intel-17.0.0-openmpi-2.0.1/install.log>` inc. install-time test results;
* :download:`Modulefile</sharc/software/modulefiles/libs/petsc/3.8.3/intel-17.0.0-openmpi-2.0.1>` installed as ``/usr/local/modulefiles/libs/petsc/3.8.3/intel-17.0.0-openmpi-2.0.1``

3.8.3 (GCC build)
^^^^^^^^^^^^^^^^^

* :download:`Install script</sharc/software/install_scripts/libs/petsc/3.8.3/gcc-6.2-openmpi-2.0.1/install.sh>`;
* :download:`Install log</sharc/software/install_scripts/libs/petsc/3.8.3/gcc-6.2-openmpi-2.0.1/install.log>` inc. install-time test results;
* :download:`Modulefile</sharc/software/modulefiles/libs/petsc/3.8.3/gcc-6.2-openmpi-2.0.1>` installed as ``/usr/local/modulefiles/libs/petsc/3.8.3/gcc-6.2-openmpi-2.0.1``
