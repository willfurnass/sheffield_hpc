OpenMPI (gcc version)
=====================

.. sidebar:: OpenMPI (gcc version)

   :Latest Version: 2.0.1
   :Dependancies: gcc
   :URL: http://www.open-mpi.org/

The Open MPI Project is an open source Message Passing Interface implementation that is developed and maintained by a consortium of academic, research, and industry partners. Open MPI is therefore able to combine the expertise, technologies, and resources from all across the High Performance Computing community in order to build the best MPI library available. Open MPI offers advantages for system and software vendors, application developers and computer science researchers.

Versions
--------

You can load a specific version using ::

   module load mpi/openmpi/2.0.1/gcc-6.2
   module load mpi/openmpi/2.0.1/gcc-5.4
   module load mpi/openmpi/2.0.1/gcc-4.9.4
   module load mpi/openmpi/1.10.4/gcc-6.2
   module load mpi/openmpi/1.10.4/gcc-4.9.4

See `here <https://mail-archive.com/announce@lists.open-mpi.org/msg00085.html>`__ for a brief guide to the new features in OpenMPI 2.x and `here <https://raw.githubusercontent.com/open-mpi/ompi/v2.x/NEWS>`__ for a detailed view of the changes between OpenMPI versions.

**CUDA**: Note that if you are using :ref:`CUDA <cuda_sharc>` with OpenMPI then you currently need to use a version of CUDA built with GCC < 5.0.
**C++ bindings** If you are using the C++ bindings then you should use OpenMPI 1.10.4 as the bindings have been deprecated in OpenMPI 2.0.1.

Examples
--------

Example programs are available in the ``$MPI_HOME/examples/`` directory.  

To compile and run these programs: copy that directory to your home directory, start an interactive MPI-aware session on a worker node, activate the version of OpenMPI you wish to use, compile the examples then run them.

In more detail ::

    # Connect to ShARC
    ssh user@sharc  

    # Start an interactive session from which we can run MPI processes using a core on each of four nodes
    qrsh -pe mpi 4

    # Load an MPI implementation
    module load mpi/openmpi/2.0.1/gcc-6.2

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

**Version 2.0.1, gcc 6.2**

1. Enable :ref:`GCC <gcc_sharc>` 6.2.0.
2. Download, compile and install OpenMPI 2.0.1 using the :download:`this script </sharc/software/install_scripts/mpi/openmpi/2.0.1/gcc-6.2/install.sh>`.
3. Install :download:`this modulefile </sharc/software/modulefiles/mpi/openmpi/2.0.1/gcc-6.2>` as ``/usr/local/modulefiles/mpi/openmpi/2.0.1/gcc-6.2``

**Version 2.0.1, gcc 5.4**

1. Enable :ref:`GCC <gcc_sharc>` 5.4.0
2. Download, compile and install OpenMPI 2.0.1 using the :download:`this script </sharc/software/install_scripts/mpi/openmpi/2.0.1/gcc-5.4/install.sh>`.
3. Install :download:`this modulefile </sharc/software/modulefiles/mpi/openmpi/2.0.1/gcc-5.4>` as ``/usr/local/modulefiles/mpi/openmpi/2.0.1/gcc-5.4``

**Version 2.0.1, gcc 4.9.4**

1. Download, compile and install OpenMPI 2.0.1 using the :download:`this script </sharc/software/install_scripts/mpi/openmpi/2.0.1/gcc-4.9.4/install.sh>`.
2. Install :download:`this modulefile </sharc/software/modulefiles/mpi/openmpi/2.0.1/gcc-4.9.4>` as ``/usr/local/modulefiles/mpi/openmpi/2.0.1/gcc-4.9.4``

**Version 1.10.4, gcc 6.2**

1. Enable :ref:`GCC <gcc_sharc>` 6.2.0.
2. Download, compile and install OpenMPI 1.10.4 using :download:`this script </sharc/software/install_scripts/mpi/openmpi/1.10.4/gcc-6.2/install.sh>`.
3. Install :download:`this modulefile </sharc/software/modulefiles/mpi/openmpi/1.10.4/gcc-6.2>` as ``/usr/local/modulefiles/mpi/openmpi/1.10.4/gcc-6.2``

	
	

	
	
 

