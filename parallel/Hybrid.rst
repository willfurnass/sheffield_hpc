.. _parallel_hybrid:

Hybrid SMP / MPI
================

Overview
--------

If you want to distribute work between a large number of cores (e.g. >16) on ShARC or Iceberg 
then the most common approach is MPI, 
as :ref:`described here <parallel_mpi>`.
The default MPI setup on ShARC and Iceberg does not ensure that your MPI job is assigned the same number of CPU cores per node.
However in certain cases you may want a symmetric distribution of of cores: 
you may want eight MPI processes, 
each of which internally parallelises work using four :ref:`OpenMP threads <parallel_smp>`, 
which would run optimally if each process were assigned four cores on a node 
(possibly running multiple processes on a node).

Usage on Iceberg
----------------

Support for Hybrid SMP/MPI is in the preliminary stages.
Here is an example job submission for an executable that combines SMP (specifically OpenMP) with MPI: ::

  #$ -pe openmpi-hybrid-4 8
  #$ -l rmem=2G
  #$ -l mem=6G
  module load mpi/intel/openmpi
  export OMP_NUM_THREADS=4
  mpirun -bynode -np 2 -cpus-per-proc 4 [executable + options]

There would be 2 MPI processes running, each with 4x (6GB mem, 2GB rmem) shared between the four threads per MPI process, for a total of 8 x (6GB mem, 2GB rmem) for the entire job.
When we tried this, we got warnings saying that the `-cpus-per-proc` was getting deprecated.  A quick google suggests that ::

  mpirun -np 2 --map-by node:pe=1 [executable + options]

would be the appropriate replacement.

Here, ``-pe openmpi-hybrid-4 8`` means 'request eight *slots* (CPU cores) from the scheduler using the ``openmpi-hybrid-4`` *Parallel Environment*'.  
This Parallel Environment ensure that the slots are assigned to nodes in groups of four.
with this example, the job could be assigned four cores on each of two nodes or eight cores on one node (given that all nodes have at least eight cores)

To see other available hybrid SMP/MPI Parallel Environments run ``qconf -spl | grep hybrid`` e.g.: ::

    openmpi-hybrid-12
    openmpi-hybrid-16
    openmpi-hybrid-4
    openmpi-hybrid-6
    openmpi-hybrid-8

Note that the number of cores you request should be wholely divisible by the number in the name of the Parallel Environment you use.

Usage on ShARC
--------------

This is similar to Iceberg but at present only one Parallel Environnment is defined: ``mpi-smp-16``.
This particular Parallel Environment can be used to request exclusive access to multiple nodes 
(as all standard nodes in ShARC have 16 cores) e.g. 
including the following in your job submission script would result in 
your job being allocated three nodes: ::

  #$ -pe mpi-smp-16 48
    
More Parallel Environments for Hybrid SMP/MPI can be added on request.
