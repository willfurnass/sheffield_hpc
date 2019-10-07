.. _impi_bessemer:

Intel MPI
=========

.. sidebar:: Intel MPI

   :URL: https://software.intel.com/en-us/mpi-library

"Intel MPI Library is a multifabric message-passing library
that implements the open-source MPICH specification.
Use the library to create, maintain, and test advanced, complex applications that
perform better on HPC clusters based on Intel® processors."

Versions
--------

You can load a specific version using: ::

    module load impi/2018.4.274-iccifort-2019.1.144-GCC-8.2.0-2.31.1 

which explicitly loads versions of icc, ifort (and GCC).

Example
-------

Consider the following source code (hello.c):

.. code-block:: c

    #include <mpi.h>
    #include <stdio.h>

    int main(int argc, char** argv) {
        // Initialize the MPI environment
        MPI_Init(NULL, NULL);

        // Get the number of processes
        int world_size;
        MPI_Comm_size(MPI_COMM_WORLD, &world_size);

        // Get the rank of the process
        int world_rank;
        MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

        // Get the name of the processor
        char processor_name[MPI_MAX_PROCESSOR_NAME];
        int name_len;
        MPI_Get_processor_name(processor_name, &name_len);

        // Print off a hello world message
        printf("Hello world from processor %s, rank %d out of %d processors\n",
               processor_name, world_rank, world_size);

        // Finalize the MPI environment.
        MPI_Finalize();
    }

MPI_COMM_WORLD (which is constructed for us by MPI) encloses all of the processes in the job, so this call should return the amount of processes that were requested for the job.

Compile your source code by using on of the following commands: ::

    mpic++ hello.cpp -o file
    mpicxx hello.cpp -o file
    mpicc hello.c -o file
    mpiCC hello.c -o file


Interactive job submission
##########################


You can run your job interactively: ::

    srun file

Your output would be something like: ::

    Hello world from processor bessemer-node001.shef.ac.uk, rank 0 out of 1 processors


This is an expected behaviour since we did not specify the number of CPU cores when requesting our interactive session.
You can request an interactive node with multiple cores (4 in this example) by using the command: ::

    srun --ntasks=4 --pty bash -i

Please note that requesting multiple cores in an interactive node depends on the availability. During peak times, it is unlikely that you can successfully request a large number of cpu cores interactively.  Therefore, it may be a better approach to submit your job non-interactively. 


Non-interactive job submission
##############################

Write a shell script (minimal example) We name the script as ‘test.sh’: ::


    #!/bin/bash
    #SBATCH --nodes=1
    #SBATCH --ntasks-per-node=40

    module use /usr/local/modulefiles/staging/eb/all/
    module load OpenMPI/3.1.3-GCC-8.2.0-2.31.1

    srun file

Maximum 40 cores can be requested.

Submit your script by using the command: ::

    sbatch test.sh

Your output would be something like: ::

    Hello world from processor bessemer-node003.shef.ac.uk, rank 24 out of 40 processors
    Hello world from processor bessemer-node002.shef.ac.uk, rank 5 out of 40 processors
    ...
    Hello world from processor bessemer-node003.shef.ac.uk, rank 31 out of 40 processors
    Hello world from processor bessemer-node003.shef.ac.uk, rank 32 out of 40 processors



