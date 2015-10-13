.. _sge-intro:

Iceberg's Queue System
======================

To manage use of the iceberg cluster, there is a queue system 
(`SoGE <https://arc.liv.ac.uk/trac/SGE>`_, a derivative of the Sun Grid Engine).

The queue system works by a user requeuesting some task, either a script or an 
interactive session, be run on the cluster and then the scheduler will take
tasks from the queue based on a set of rules and priorities.


.. _sge-interactive:

Running an Interactive Shell
----------------------------

If you wish to use the cluster for interactive use, such as running applications
such as MATLAB or Ansys, or compiling software, you will need to requeuest that
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
combined together to requeuest more resources.

.. note::

    Long running jobs *should* use the batch submission system rather than 
    requeuesting an interactive session for a very long time. Doing this will 
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

.. _sge-batch:

Submitting Jobs to the queue
----------------------------

The power of iceberg really comes from the 'batch job' queue submission process.
Using this system, you write a script which executes your job, tell the 
scheduler how many resources the task requires, then the scheduler will run it 
when the resources are available.
As the task is running, the terminal output and any errors are captured and 
saved to a disk, so that you can see the output and verify the execution of the
task.

Any task that can be executed without any user intervention while it is running 
can be submitted as a batch job to iceberg. This exculdes jobs that require a 
GUI, however, many common applications such as Ansys or MATLAB can also be 
used without their GUIs.

When you submit a batch job, you provide an executable file that will be run by
the scheduler. This is normally a script file which provides commands and
options to the program you are using. For instance, it might tell Ansys which 
files to use as input and where to save the output. Once you have a script 
file, or other executable file, you can submit it to the queue by running::

    qsh myscript.sh

you can also specify extra arguments to this, or at the start of your script, 
to give you access to more cores or memory or change the maximum execution time,
a full list of the availble options are given below.


All Scheduler Options
---------------------


====================== ========================================================
Command                Description
====================== ========================================================
-l h_rt=hh:mm:ss       Specify the total maximum execution time for the job.

-l mem=xxG             Specify the maximum amount (xx) of memory to be used. 

-N                     Job name, used to name output files and in the queue list.

-j                     Join the error and normal output into one file rather 
                       than two.

-M                     Email address to send notifications to.

-m bea                 Type of notifications to send. Can be any combination of
                       begin (b) end (e) or abort (a) i.e. `-m ea` for end and 
                       abortion messages.
-a                     Specify the earliest time for a job to start, in the
                       format MMDDhhmm. e.g. -a 01011130 will schedule the job
                       to begin no sooner than 11:30 on 1st January.
====================== ========================================================

Frequently Asked SGE Questions
------------------------------
**How do you ensure that a job starts after a specified time?**

Add the following line in your submission script ::

    #$ -a time

but replace ``time`` with a time in the format MMDDhhmm

For example, for 22nd July at 14:10, you’d do ::

    #$ -a 07221410

This won’t guarantee that it will run precisely at this time since that depends on available resources. It will, however, ensure that the job runs *after* this time. If your resource requirements aren’t too heavy, it will be pretty soon after. When I tried it, it started about 10 seconds afterwards but this will vary.
