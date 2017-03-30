
Scheduler
=========

This is a reference section for commands that allow you to interact with the scheduler.

**qhost**

`qhost` is a scheduler command that show's the status of Sun Grid Engine hosts.

**qrsh**


`qrsh` is a scheduler command that requests an interactive session on a worker node. The resulting session will **not** support graphical applications. You will usually run this command from the head node.

Examples

Request an interactive session that provides the default amount of memory resources ::

    qrsh

Request an interactive session that provides 10 Gigabytes of real and virtual memory ::

    qrsh -l rmem=10G -l mem=10G

**qrshx**

`qrshx` is a scheduler command that requests an interactive session on a worker node. The resulting session will support graphical applications. You will usually run this command from the head node.

Examples

Request an interactive X-Windows session that provides the default amount of memory resources and launch the `gedit` text editor ::

    qrshx
    gedit

Request an interactive X-Windows session that provides 10 Gigabytes of real and virtual memory and launch the latest version of MATLAB ::

    qrshx -l mem=10G -l rmem=10G
    module load apps/matlab
    matlab

Request an interactive X-Windows session that provides 10 Gigabytes of real and virtual memory and 4 CPU cores ::

    qrshx -l rmem=10G -l mem=10G -pe openmp 4

Sysadmin notes

qrshx is a Sheffield-developed modification to the standard set of scheduler commands. It is at `/usr/local/bin/qrshx` and contains the following ::

  #!/bin/sh
  exec qrsh -v DISPLAY -pty y "$@" bash

**qsh**

`qsh` is a scheduler command that requests an interactive X-windows session on a worker node. The resulting terminal is not user-friendly and we recommend that you use our `qrshx` command instead.

Examples

Request an interactive X-Windows session that provides the default amount of memory resources ::

    qsh

Request an interactive X-Windows session that provides 10 Gigabytes of real and virtual memory ::

    qsh -l rmem=10G -l mem=10G

Request an interactive X-Windows session that provides 10 Gigabytes of real and virtual memory and 4 CPU cores ::

   qsh -l rmem=10G -l mem=10G -pe openmp 4


**qstat**

`qstat` is a scheduler command that displays the status of the queues.

Examples

Display all jobs queued on the system ::

    qstat

Display all jobs queued by the username foo1bar ::

    qstat -u foo1bar

Display all jobs in the openmp parallel environment ::

    stat -pe openmp

Display all jobs in the queue named foobar ::

    qstat -q foobar.q

**qsub**

`qsub` is a scheduler command that submits a batch job to the system.

Examples

Submit a batch job called myjob.sh to the system ::

    qsub myjob.sh

**qtop**

`qtop` is a scheduler command that provides a summary of all processes running on the cluster for a given user.


Examples

`qtop` is only available on the worker nodes. As such, you need to start an interactive session on a worker node using `qrsh` or `qrshx` in order to use it.

To give a summary of all of your currently running jobs ::

    qtop

    Summary for job 256127

        HOST      VIRTUAL-MEM       RSS-MEM    %CPU    %MEM    CPUTIME+   COMMAND
    testnode03      106.22 MB        1.79 MB     0.0     0.0     00:00:00 bash
    testnode03      105.62 MB        1.27 MB     0.0     0.0     00:00:00 qtop
    testnode03       57.86 MB        3.30 MB     0.0     0.0     00:00:00 ssh
                    ---------       --------
        TOTAL:        0.26 GB        0.01 GB


Scheduler Options
-----------------

====================== ========================================================
Command                Description
====================== ========================================================
-l h_rt=hh:mm:ss       Specify the total maximum execution time for the job.

-l mem=xxG             Specify the maximum amount (xx) of memory to be used.

-l hostname=           Target a node by name. Not recommended for normal use.

-l arch=               Target a processor architecture. Options on Iceberg include
                       `intel-e5-2650v2` and `intel-x5650`

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