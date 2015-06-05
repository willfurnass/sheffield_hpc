R
=

.. sidebar:: R
   
   :Support Level: bronze
   :Dependancies: BLAS
   :URL: http://www.r-project.org/ 
   :Documentation: http://www.r-project.org/  

R is a statistical computing language.

Interactive Usage
-----------------
After connecting to iceberg (see :ref:`ssh`),  start an interactive sesssion with the :code:`qsh` command.

The lastest version of R can be loaded with ::

        module load apps/R

Alternatively, you can load a specific version of R using one of the following ::

        module load apps/R/3.1.2
        module load apps/R/3.2.0

R can then be run with ::

        $ R

Serial (one CPU) Batch usage
----------------------------
Here, we assume that you wish to run the program :code:`my_code.R` on the system.

First, you need to write a submission file. We assume you'll call this :code:`my_job.sge` ::

	#!/bin/bash
	#$ -S /bin/bash
	#$ -cwd               # Run job from current directory
	
        module load apps/R/3.2.0 # Load version 3.2.0 of R
    
	R CMD BATCH my_code.R my_code.R.o$JOB_ID

Note that R must be called with both the :code:`CMD` and :code:`BATCH` options which tell it to run an R program, in this case :code:`my_code.R`. If you do not do this, R will attempt to open an interactive prompt.

The final argument, :code:`my_test.R.o$JOBID`, tells R to send output to a file with this name. Since :code:`$JOBID` will always be unique, this ensures that all of your output files are unique. Without this argument R sends all output to a file called :code:`my_code.Rout`.

Ensuring that :code:`my_code.R` and :code:`my_job.sge` are both in your current working directory, submit your job to the batch system ::

	qsub my_job.sge

Replace :code:`my_job.sge` with the name of your submission script.

Graphical output
----------------
By default, graphical output from batch jobs is sent to a file called :code:`Rplots.pdf`

Installation Notes
------------------

R was compiled from source using the following commands::

        $ module load libs/gcc/lapack
        $ module load libs/gcc/blas
        $ ./configure --prefix /usr/local/packages6/R/3.1.2 --use-blas --use-lapack
        $ make -j 6
        $ make install

It seems like it did not pick up the lapack library, but it did pick up the BLAS lib ok.
