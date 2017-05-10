.. _bclr:

Berkeley Lab Checkpoint/Restart (BLCR)
======================================

.. sidebar:: BLCR

   :Latest version: 0.8.5
   :URL: http://crd.lbl.gov/departments/computer-science/CLaSS/research/BLCR/

From the `BLCR User's Guide <https://upc-bugs.lbl.gov/blcr/doc/html/BLCR_Users_Guide.html>`_: 

Checkpoint/Restart allows you to **save one or more processes to a file and later restart them from that file**. There are three main uses for this:
 
#. Scheduling: Checkpointing a program allows a program to be **safely stopped at any point in its execution**, so that some other program can run in its place. The original program can then be run again later.
#. Process Migration: If a compute node appears to be likely to crash, or there is some other reason for shutting it down (routine maintenance, hardware upgrade, etc.), checkpoint/restart allows any processes running on it to be moved to a different node (or saved until the original node is available again).
#. Failure recovery: A **long-running program can be checkpointed periodically, so that if it crashes due to hardware, system software, or some other non-deterministic cause, it can be restarted from a point in its execution more recent that starting from the beginning**.
 
Berkeley Lab Checkpoint/Restart (BLCR) provides checkpoint/restart on Linux systems. BLCR can be used either with processes on a single computer, or on parallel jobs (such as MPI applications) which may be running across multiple machines on a cluster of Linux nodes.
 
Note: Checkpointing parallel jobs requires a library which has integrated BLCR support. At the present time, many MPI implementations are known to support checkpoint/restart with BLCR. Consult the corresponding BLCR `entry <https://upc-bugs.lbl.gov/blcr/doc/html/FAQ.html#mpi>`_ for the current list. 

Usage
-----
BCLR was installed before we started using the module system. As such, it is currently always available to worker nodes. It should be considered experimental.

The checkpointing is performed at the kernel level, so any batch code should be checkpointable without modification (it may not work with our MPI environment though, although it should cope with SMP codes).

To run a code, use ::

    cr_run ./executable

To checkpoint a process with process id PID ::

    cr_checkpoint -f checkpoint.file PID

Use the ``--term`` flag if you want to checkpoint and kill the process

To restart the process from a checkpoint ::

    cr_restart checkpoint.file

Using with SGE in batch
-----------------------
A checkpoint environment has been setup called ``blcr``.  An example of a checkpointing job would look something like ::

    #!/bin/bash
    #$ -l h_rt=168:00:00
    #$ -c sx
    #$ -ckpt blcr
    cr_run ./executable >> output.file

The ``-c sx`` options tells the scheduler to checkpoint if the queue is suspended, or if the execution daemon is killed. You can also specify checkpoints to occur after a given time period.

A checkpoint file will be produced before the job is terminated.  This file will be called ``checkpoint.[jobid].[pid]``.  This file will contain the complete in-memory state of your program at the time that it terminates, so make sure that you have enough disk space to save the file.

To resume a checkpointed job submit a job file which looks like ::

  #!/bin/bash
  #$ -l h_rt=168:00:00
  #$ -c sx
  #$ -ckpt blcr
  cr_restart [name of checkpoint file]

Installation notes
------------------
Installation notes are not available.
