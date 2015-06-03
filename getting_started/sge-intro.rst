.. _sge-intro:

iceberg's Que System
====================

To manage use of the iceberg cluster, there is a que system 
(`SoGE <https://arc.liv.ac.uk/trac/SGE>`_, a derivative of the Sun Grid Engine).

The que system works by a user requesting some task, either a script or an 
interactive session, be run on the cluster and then the scheduler will take
tasks from the que based on a set of rules and priorities.


.. _sge-interactive:

Running an Interactive Shell
############################

If you wish to use the cluster for interactive use, such as running applications
such as MATLAB or Ansys, or compiling software, you will need to request that
the scheduler gives you an interactive session. For an introduction to this see
:ref:`getting-started`.

There are two commands which give you an interactive shell::

    [te1st@iceberg-login2 ~]$ qsh

and::
    
    [te1st@iceberg-login2 ~]$ qrsh

qsh will open a separate xterm terminal and supports running graphical 
applications. qrsh gives a shell running on a worker node inside the currently 
open terminal, it does not support graphical applications because it has no 
X server forwarding configured.

You can configure the resources available to the interactive session by 
specifying them as command line options to the qsh or qrsh commands.
For example to run a qsh session with access to 16 GB of virtual RAM::


    [te1st@iceberg-login2 ~]$ qsh -l mem=16G

or a session with access to 8 cores::


    [te1st@iceberg-login2 ~]$ qsh -pe openmp 8

A table of :ref:`sge-interactive-options` is given below, any of these can be 
combined together to request more resources.

.. note::

    Long running jobs *should* use the batch submission system rather than 
    requesting an interactive session for a very long time. Doing this will 
    lead to better cluster performance for all users.


.. _sge-interactive-options:

Common Interactive Job Options
``````````````````````````````

====================== ========================================================
Command                Description
====================== ========================================================
-l h_rt=hh:mm:ss       Specify the total maximum execution time for the job.

-l mem=xxG             Specify the maximum amount (xx) of memory to be used 
                       (per process or core). 

-pe <env> <nn>         Specify a parallel environment and number of processors. 

-pe openmp <nn>        The openmp parallel environment provides multiple threads
                       on one node. It is the most commonly used in an 
                       interactive session. <nn> specifies the max number of 
                       threads.

-pe openmpi-ib <nn>    The openmpi parallel environment is for distributed 
                       memory computing, it enables multiple *processes* to 
                       run independently. It is more commonly used in batch 
                       mode.
====================== ========================================================

Submitting Jobs to the Que
##########################

The power of iceberg really comes from the 'batch job' que submission process.
Using this system, you write a script which executes your job, tell the 
scheduler how many resources the task requires, then the scheduler will run it 
when the resources are available.

As the task is running, the terminal output and any errors are captured and 
saved to a disk, so that you can see the output and verify the execution of the
task.


All Scheduler Options
#####################


====================== ========================================================
Command                Description
====================== ========================================================
-l h_rt=hh:mm:ss       Specify the total maximum execution time for the job.

-l mem=xxG             Specify the maximum amount (xx) of memory to be used. 

-N                     Job name, used to name output files and in the que list.

-j                     Join the error and normal output into one file rather 
                       than two.

-M                     Email address to send notifications to.

-m bea                 Type of notifications to send. Can be any combination of
                       begin (b) end (e) or abort (a) i.e. `-m ea` for end and 
                       abortion messages.
====================== ========================================================

