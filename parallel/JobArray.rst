.. _parallel_jobarray:

SGE Job Arrays
==============

The simplest way of exploiting parallelism on the clusters is to use **Job Arrays**. A Job Array is a set of batch jobs run from a single job script. For example: ::

  #!/bin/bash
  #
  #$ -t 1-100
  #
  echo "Task id is $SGE_TASK_ID"

  ./myprog.exe $SGE_TASK_ID > output.$SGE_TASK_ID

The above script will submit 100 tasks to the system at once.
The difference between each of these 100 jobs is the value of the environment variable ``$SGE_TASK_ID`` which will range from 1 to 100, determined by the line ``#$ -t 1-100``.
In this example, the program ``myprog.exe`` will be run 100 times with differing input values of ``$SGE_TASK_ID``. 100 output files will be created with filenames ``output.1``, ``output.2`` and so on.

Job arrays are particularly useful for `Embarrassingly Parallel <https://en.wikipedia.org/wiki/Embarrassingly_parallel>`_ problems such as Monte Carlo Simulations (Where ``$SGE_TASK_ID`` might correspond to random number seed), or batch file processing (where ``$SGE_TASK_ID`` might refer to a file in a list of files to be processed).

Examples
--------
* `MATLAB SGE Job Array example <https://github.com/mikecroucher/HPC_Examples/tree/master/languages/MATLAB/SGE_array>`_

Email notifications
-------------------

If you :ref:`enable email notifications <sge-queue>` in your batch job submission script then you will receive emails for *every* job in your job array.  This helps you determine if any jobs in the array failed but doesn't help you determine if the entire job array has finished running.  Here's a sensible approach to managing email notifications:

1. Edit your job array submission script so you are *only* notified of aborted (``a``) jobs i.e. ::

        #$ -M me@sheffield.ac.uk
        #$ -m a

2. Then submit your job array like so: ::

        [te1st@sharc-login1 ~]$ qsub my_job_array.sge
        Your job-array 2035587.1-3:1 ("my_job_array.sge") has been submitted

3. Next, submit a very simple job that will only run when job array ``2035587`` has completed and emails you when it finishes: ::

        [te1st@sharc-login1 ~]$ qsub -o /dev/null -e /dev/null -M me@sheffield.ac.uk -m ea -b y -l h_rt=00:00:15 -hold_jid 2035587 -N 'Job_array_finished' true
        Your job 2035588 ("Job_array_finished") has been submitted

You will therefore receive:

* An email for every failed job in the array
* An email shortly after the entire job array finishes
