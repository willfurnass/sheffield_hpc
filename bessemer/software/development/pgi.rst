.. _`PGI Compilers_bessemer`:

PGI Compilers
=============

The PGI Compiler suite offers C, C++ and Fortran Compilers.
For full details of the features of this compiler suite
see `PGI's website <http://www.pgroup.com/products/pgiworkstation.htm>`_.

Making the PGI Compilers available
----------------------------------

After connecting to the Bessemer cluster, start an interactive session: ::

   srun --pty bash -i

then activate a specific version of the compiler suite using: ::

   module load PGI/19.1-GCC-8.2.0-2.31.1

Once you've loaded the module, you can check the version with: ::

   pgcc --version


Compilation examples
--------------------


**C**

To compile a C hello world example into an executable called ``hello`` using the PGI C compiler ::

    pgcc hello.c -o hello

**C++**

To compile a C++ hello world example into an executable called ``hello`` using the PGI C++ compiler ::

    pgc++ hello.cpp -o hello

**Fortran**

To compile a Fortran hello world example into an executable called ``hello`` using the PGI Fortran compiler ::

    pgf90 hello.f90 -o hello


Interactive Usage
-----------------

Run the compiled program by using the command `./file`.


Non-Interactive Usage
---------------------

Write a shell script (minimal example) We name the script as ‘test.sh’: ::

    #!/bin/bash
    srun file

Submit your script by using the command: ::

    sbatch test.sh

Note the job submission number,  for example: Submitted batch job 1211

Check your output file:  ::

    cat slurm-1211.out

Additional options:

Name your submission: ::
    #SBATCH --comment=test_job

Specify nodes and tasks: ::
    #SBATCH --nodes=1
    #SBATCH --ntasks-per-node=16

Memory allocation: ::
    #SBATCH --mem=16000

Specify the output file name: ::
    #SBATCH --output=output.%j.test.out

Request time: ::
    #SBATCH --time=00:30:00

Email notification: ::
    #SBATCH --mail-user=username@sheffield.ac.uk

