.. _sched_delete:

Cancelling jobs
===============

Deleting queued/running jobs can be useful/necessary if you realise that:

- there is an issue with the software and/or data;
- you have requested too much in the way of resources and the job will take a long time to be scheduled;
- you have requested way too much in the way of resources and the job will never run;
- you have requested too little in the way of resources and the job will probably be killed.

.. contents::
   :local:

On SGE
------

You can delete a job by first finding the job ID using ``qstat`` (see :ref:`queued_running`) then running: ::

    qdel -j myjobid

Alternatively you can delete all of your queued/running jobs using: ::

    qdel -a

``qdel`` has many more options; to read up on these log in to the cluster and run: ::

    man qdel


On Slurm
--------

You can delete a job by first finding the job ID using ``sacct -u $USER`` (see :ref:`queued_running`) then running: ::

    scancel myjobid

Alternatively you can delete all of your queued/running jobs using: ::

    scancel -u $USER

``scancel`` has many more options; to read up on these log in to the cluster and run: ::

    man scancel

