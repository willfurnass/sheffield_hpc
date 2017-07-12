.. _mvapich2_intel_sharc:

MVAPICH2 (Intel version)
=======================

.. sidebar:: MVAPICH2 (Intel version)

   :Latest Version: 2.3a
   :Dependancies: Intel compilers
   :URL: http://mvapich.cse.ohio-state.edu/

MVAPICH, also known as MVAPICH2, is, like OpenMPI, an implementation of the MPI standard for passing messages between machines in parallel computing environments.

Versions
--------

You can load MVAPICH2 using ::

   module load mpi/mvapich/2.3a/intel-17.0.0

This version was build with the Intel compiler suite version 17.0.0.

See `here <http://mvapich.cse.ohio-state.edu/overview/>`__ for an overview of the features offered per MVAPICH2 release.

Examples - ADD MVAPICH2 EXAMPLES OR REMOVE THIS SECTION
-------------------------------------------------------

Example programs are available in the ``/usr/local/packages/mpi/openmpi/XXX/intel-17.0.0/examples/`` directory (where ``XXX`` is the version of OpenMPI you are using e.g. ``2.0.1``).  

To compile and run these programs: copy that directory to your home directory, start an interactive MPI-aware session on a worker node, activate the version of OpenMPI you wish to use, compile the examples then run them.

In more detail ::

    # Connect to ShARC
    ssh user@sharc  

    # Start an interactive session from which we can run MPI processes using a core on each of four nodes
    qrsh -pe mpi 4

    # Load an MPI implementation
    module load mpi/openmpi/2.0.1/intel-17.0.0

    # Copy the examples to your home directory
    cp -r $MPI_HOME/examples ~/openmpi_2.0.1_examples

    # Compile all programs in the examples directory
    cd ~/openmpi_2.0.1_examples
    make

    # Once compiled, run an example program on all (or a subset) of your MPI nodes using the mpirun utility
    mpirun -np 4 hello_c
    

    Hello, world, I am 0 of 4, (Open MPI v2.0.1, package: Open MPI user@sharc-node002.shef.ac.uk Distribution, ident: 2.0.1, repo rev: v2.0.0-257-gee86e07, Sep 02, 2016, 141)
    Hello, world, I am 1 of 4, (Open MPI v2.0.1, package: Open MPI user@sharc-node002.shef.ac.uk Distribution, ident: 2.0.1, repo rev: v2.0.0-257-gee86e07, Sep 02, 2016, 141)
    Hello, world, I am 2 of 4, (Open MPI v2.0.1, package: Open MPI user@sharc-node002.shef.ac.uk Distribution, ident: 2.0.1, repo rev: v2.0.0-257-gee86e07, Sep 02, 2016, 141)
    Hello, world, I am 3 of 4, (Open MPI v2.0.1, package: Open MPI user@sharc-node002.shef.ac.uk Distribution, ident: 2.0.1, repo rev: v2.0.0-257-gee86e07, Sep 02, 2016, 141)


Installation notes
------------------

These are primarily for administrators of the system.

Version 2.3a, Intel 17.0.0 compiler
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

#. Ensure :ref:`Intel compilers 17.0.0 <sharc-intel-compilers>` are installed and licensed.
#. Download, compile and install MVAPICH2 2.3a using :download:`this script </sharc/software/install_scripts/mpi/mvapich2/2.3a/intel-17.0.0/install.sh>`.
#. Install :download:`this modulefile </sharc/software/modulefiles/mpi/mvapich2/2.3a/intel-17.0.0>` as ``/usr/local/modulefiles/mpi/mvapich2/2.3a/intel-17.0.0``
#. Test by running some <Examples>_. **ADD EXAMPLES OR REMOVE THIS LINE**
