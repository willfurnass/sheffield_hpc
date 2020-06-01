.. _sched:

Doing work using interactive sessions / batch jobs
==================================================

If you want to do computational work on a cluster then you need to either:

- Start an **interactive session** on one (or more) worker node(s) or
- Submit a **batch job** that will run on one (or more) worker node(s)

An **interactive session** (typically) allows you to interactively run programs from the command-line.
This is useful for running graphical programs (e.g. commercial engineering applications),
for debugging (where interactivity is required) and for shorter, computationally-intensive tasks (e.g. compiling a program).
When requesting an interactive session you can state the resources (CPU cores, memory, GPUs) you want for the session (including the duration)
and if the resources are available the session will start straight away, typically with exclusive access to the requested resources.

A **batch job** is started using a script that defines the things required to run a computational task at some point in the future:
in the script you specify the resources the job requires and also state the commands that should be run
(e.g. to ensure certain versions of software are available or to start a particular Python program).
After creating a batch job script, you submit it to software known as the **job scheduler**.
The job, along with other jobs submitted by you and other users, will then queue whilst waiting for
sufficient resources to become free to run your job.  Your priority in the queue depends on various factors including how 'big' your job is and how much you have been using the cluster recently.
When your job reaches the front of a queue and starts running the job scheduler monitors resource utilisation (and checks the initial resource request is not exceeded).
When the job finishes the scheduler can report to the user whether the job finished successfully and what resources were used whilst it was active.

**Users are encouraged encouraged to perform work using batch jobs where possible** as:

- Non-trivial resource requests are difficult to satisfy on-demand
- It leads to better cluster utilisation overall (fewer jobs are waiting on a user to type something)

Interactive and batch jobs at Sheffield are managed using one of two job scheduling systems:
the `Son of Grid Engine <https://arc.liv.ac.uk/trac/SGE>`_ **job scheduling software**
(often referred to as **SGE**, as it is one of several derivatives of `Sun Grid Engine <https://en.wikipedia.org/wiki/Oracle_Grid_Engine>`_)
and `Slurm <https://slurm.schedmd.com>`_ .

+----------+-----------+
| Cluster  | Scheduler |
+==========+===========+
| Bessemer | Slurm     |
+----------+-----------+
| ShARC    | SGE       |
+----------+-----------+
| Iceberg  | SGE       |
+----------+-----------+

Both schedulers work as follows: a user requests that a *job* (task), either a script or an
interactive session, be run on the cluster and then the scheduler will take jobs from
the queue based on a set of rules and priorities.

.. toctree::
    :maxdepth: 1
    :glob:

    ./interactive
    ./batch
    ./options
    ./monitor
    ./cancel
    ./history
    ./node_info
    ./programmatic

.. todo:
   add section on altering queued jobs
