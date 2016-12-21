.. _gulp_sharc:

GULP
====

.. sidebar:: GULP

   :Latest Version: 4.4
   :URL: https://nanochemistry.curtin.edu.au/gulp/

Overview
--------

GULP is a program for performing a variety of types of simulation on materials
using boundary conditions of 0-D (molecules and clusters), 1-D (polymers), 2-D
(surfaces, slabs and grain boundaries), or 3-D (periodic solids). The focus of
the code is on analytical solutions, through the use of lattice dynamics, where
possible, rather than on molecular dynamics. A variety of force fields can be
used within GULP spanning the shell model for ionic materials, molecular
mechanics for organic systems, the embedded atom model for metals and the
reactive REBO potential for hydrocarbons. Analytic derivatives are included up
to at least second order for most force fields, and to third order for many.

GULP can either be run:

- on one core (*serial* build) or
- on multiple cores, potentially on different machines (*MPI* build).

Licensing
---------

The authors of GULP have given the University of Sheffield permission to install the software on our clusters.  
The license restrictions are as follows:

 * The program is not to be distributed to anyone else without the express permission of the author.
 * The program is **not to be used for commercial research**. For any commercial use of the program a license must be obtained from Accelrys Inc, including contract research.
 * The program is supplied on an "as is" basis with no implied guarantee or support.

To ensure that the software is not used for commercial purposes access to GULP is restricted to members of the ``gulp`` (UNIX) group.  
To be added to this group, you will need to contact ``research-it@sheffield.ac.uk`` and state the nature of your work.

Interactive Usage
-----------------

The serial build of GULP should be used for interactive usage. 
After connecting to ShARC (see :ref:`ssh`),  
start an interactive session with the ``qrsh`` or ``qsh`` command. 
Make the serial version of GULP available using the following: ::

        module load apps/gulp/4.4/intel-17.0.0

You can then run GULP and show all output on the screen using: ::

        gulp < inputfile

You can also simultaneously save this output in a file (here called ``outputfile``): ::

        gulp < inputfile | tee outputfile

If you get: ::

        -bash: gulp: command not found

it is probably because you are not a member of the ``gulp`` group. See Licensing_.

Serial batch jobs
-----------------

Create a batch job submission script like the following: ::

        #!/bin/bash
        #$ -N my_gulp_serial_job
        #$ -j y
        #$ -m bea
        #$ -M YOUREMAIL@sheffield.ac.uk

        module load apps/gulp/4.4/intel-17.0.0

        gulp < infile > outfile

Then submit this job to the scheduler: ::

        qsub my_batch_job_script.sge

Multi-node batch jobs
---------------------

Create a batch job submission script like the following: ::

        #!/bin/bash
        #$ -N my_gulp_mpi_16_slot_job
        #$ -pe mpi 16
        #$ -j y
        #$ -m bea
        #$ -M YOUREMAIL@sheffield.ac.uk

        module load apps/gulp/4.4/intel-17.0.0-openmpi-2.0.1

        mpirun -np $NSLOTS gulp < infile > outfile

Then submit this job to the scheduler: ::

        qsub my_batch_job_script.sge

Examples
--------

The software comes with several example input files (plus example output files for verifying that the software is working correctly).

To run them, first activate either serial or the MPI build of GULP (using ``module load``) then copy the examples to a writable location e.g.: ::
        
        cp $GULP_EXAMPLES /data/$USER/gulp_4_4_examples
        cd /data/$USER/gulp_4_4_examples

Next, create a batch job submission script like the following (for serial testing): ::

        #!/bin/bash
        #$ -N gulp_ex_serial
        #$ -j y
        #$ -m bea
        #$ -M YOUREMAIL@sheffield.ac.uk

        module load apps/gulp/4.4/intel-17.0.0

        for infile in example*.gin; do
            echo "Running ${infile}:"
            outfile=${infile/gin/got}
            time (gulp < $infile > $outfile)
        done
        chmod +x ./diff.sh
        ./diff.sh

or like the following (for MPI testing using 16 cores): ::

        #!/bin/bash
        #$ -N gulp_ex_mpi_16
        #$ -pe mpi 16
        #$ -j y
        #$ -m bea
        #$ -M YOUREMAIL@sheffield.ac.uk

        module load apps/gulp/4.4/intel-17.0.0-openmpi-2.0.1

        for infile in example*.gin; do
            echo "Running ${infile}:"
            outfile=${infile/gin/got}
            time (mpirun -np $NSLOTS gulp < $infile > $outfile)
        done
        chmod +x ./diff.sh
        ./diff.sh

Finally, submit this job to the scheduler: ::

        qsub my_batch_job_script.sge

After receiving email notification that the job has finished, check in the ``gulp_4_4_examples`` directory for an output file containing:

 - the names of the examples that were run;
 - timings per example;
 - details of any errors;
 - details of any differences between the generated outputs and the sample outputs provided by the software's authors.

Documentation
-------------

See the `GULP website <https://nanochemistry.curtin.edu.au/gulp/>`_ or the files in the ``$GULP_DOCS`` directory on the cluster.

Installation Notes
------------------

These are primarily for system administrators.

Version 4.4
^^^^^^^^^^^

Serial (1 CPU core) and parallel (MPI) builds were compiled. 
Both builds were compiled with :ref:`version 17.0.0 of the Intel Fortran compiler <sharc-intel-compilers>` and the :ref:`Intel MKL 2017.1 <sharc-intel-mkl>`.
The MPI build was compiled using :ref:`OpenMPI 2.0.1 <openmpi_intel_sharc>`.

Both builds were compiled and installed using :download:`this script </sharc/software/install_scripts/apps/gulp/4.4/intel-17.0.0/install.sh>` plus

* :download:`this serial build configuration</sharc/software/install_scripts/apps/gulp/4.4/intel-17.0.0/getmachine_serial>`
* :download:`this MPI build configuration</sharc/software/install_scripts/apps/gulp/4.4/intel-17.0.0/getmachine_mpi>`

In addition:

* :download:`The serial build modulefile </sharc/software/modulefiles/apps/gulp/4.4/intel-17.0.0>` was installed as ``/usr/local/modulefiles/apps/gulp/4.4/intel-17.0.0``
* :download:`The parallel build modulefile </sharc/software/modulefiles/apps/gulp/4.4/intel-17.0.0-openmpi-2.0.1>` was installed as ``/usr/local/modulefiles/apps/gulp/4.4/intel-17.0.0-openmpi-2.0.1``

Both versions were tested using the examples provided by GULP's authors.  Timings per example are saved here:

* **serial**: ``/usr/local/packages/apps/gulp/4.4/intel-17.0.0/Examples/timings_2016-11-25.txt``
* **MPI**: ???
