.. _mvapich2_intel_sharc:

MVAPICH2 (Intel version)
=======================

.. sidebar:: MVAPICH2 (Intel version)

   :Latest Version: 2.3b
   :Dependancies: Intel compilers
   :URL: http://mvapich.cse.ohio-state.edu/

MVAPICH, also known as MVAPICH2, is, like OpenMPI, an implementation of the MPI standard for passing messages between machines in parallel computing environments.

Versions
--------

You can load MVAPICH2 using ::

   module load mpi/mvapich/2.3b/intel-17.0.0

This version was build with the Intel compiler suite version 17.0.0.

See `here <http://mvapich.cse.ohio-state.edu/overview/>`__ for an overview of the features offered per MVAPICH2 release.

.. _mvapich2_benchmark_progs:

Running MPI programs
--------------------

A set of benchmarking programs programs, the *OSU Micro Benchmarks*, are available in the directory: ::

    ${MPI_HOME}/libexec/osu-micro-benchmarks/mpi/

Documentation and source code for these benchmarks can be found `online <http://mvapich.cse.ohio-state.edu/benchmarks/>`__.

To run one of these programs: start an interactive MPI-aware session on a worker node, load the MVAPICH2 module file, compile the examples then run them.

In more detail ::

    # Connect to ShARC
    ssh user@sharc  

    # Start an interactive session from which we can run MPI processes using a core on each of four nodes
    qrsh -pe mpi-rsh 4

    # Display the scheduler's allocation of 'slots' (effectively CPU cores) for this job
    echo $PE_HOSTFILE

    # Load an MPI implementation
    module load mpi/mvapich2/2.3b/intel-17.0.0

    # run an benchmark program on all of your MPI nodes using the mpirun utility
    mpirun hello_c
    # or to run on a subset of the nodes
    mpirun -np 2 hello_c

Installation notes
------------------

These are primarily for administrators of the system.

Version 2.3b, Intel 17.0.0 compiler
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

#. Ensure :ref:`Intel compilers 17.0.0 <sharc-intel-compilers>` are installed and licensed.
#. Download, compile and install MVAPICH2 2.3b using :download:`this script </sharc/software/install_scripts/mpi/mvapich2/2.3b/intel-17.0.0/install.sh>`.
#. Install :download:`this modulefile </sharc/software/modulefiles/mpi/mvapich2/2.3b/intel-17.0.0>` as ``/usr/local/modulefiles/mpi/mvapich2/2.3b/intel-17.0.0``
#. Test by running :ref:`OSU micro benchmarks <mvapich2_benchmark_progs>'.
