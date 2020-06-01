.. _sched_options:

Job submission options
======================

Some of the more common/useful options that can be used when
starting :ref:`interactive sessions <sched_interactive>` or :ref:`batch job <sched_batch>`.
Note that in the table below ``XX`` and ``YY`` are placeholders that should be replaced with values of your choosing.

.. list-table::
   :widths: auto
   :align: left
   :header-rows: 1

   * - SGE
     - Slurm
     - Description
   * - ``-l h_rt=hh:mm:ss``
     - ``-t [min]`` or ``-t [days-hh:mm:ss]``
     - Specify the total maximum execution time for the job.
       The upper limit is typically 96:00:00 (4 days) on ShARC
       and 168:00:00 (7 days) on Iceberg and Bessemer.
       Note that these limits may differ for specific Projects/Queues/Partitions.
       Also note that requesting less execution time may result in your job spending less time queuing.

   * - ``-pe <env> XX``
     - N/A
     - Specify a *parallel environment*
       and a number of processor cores.

   * - ``-pe smp XX``
     - ``-c XX``
     - For parallel jobs requesting XX CPU cores on a single node

   * - ``-l rmem=XXG``
     - ``--mem=XXG``
     - Specify the maximum amount of real memory to be requested
       **per requested CPU core**.
       If the real memory usage of your job exceeds this value
       multiplied by the number of cores you requested
       then your job will be killed.

   * - ``-l arch=XX``
     - N/A
     - Target a processor architecture.
       This is irrelevant if using ShARC's *public* nodes
       as all processors are the same model.
       Options on Iceberg include ``intel-e5-2650v2`` and ``intel-x5650``.

   * - ``-l gpu=XX``
     - ``--gres=gpu:XX``
     - Request XX GPUs in a single node.
       With SGE this request is *per requested CPU core*.

   * - ``-l gpu_arch=XX``
     - TBC
     - Request a particular model of GPU.

   * - ``-N XX``
     - ``--job-name=XX``
     - Job name; shown in job queue information and used by default in names of output files

   * - ``-o XX``
     - ``--output=XX``
     - File to send the normal job output (``stdout``) to.
       On Slurm this can be templated using special variables
       to insert e.g. the job ID; run ``man sbatch`` for more info.

   * - ``-e XX``
     - ``--error=XX``
     - File to send the error job output (``stderr``) to.
       On Slurm this can be templated using special variables
       to insert e.g. the job ID; run ``man sbatch`` for more info.

   * - ``-j y[es]|n[o]``
     - ``-o [filename]``
     - Join the error output and normal output in one file rather
       than send them to separate files.

   * - ``-M XX``
     - ``--mail-user=XX``
     - Email address to send notifications to.

   * - ``-m bea``
     - ``--mail-type=XX``
     - Type of notifications to send.
       For SGE can be any combination of
       begin (``b``) end (``e``) or abort (``a``) i.e.
       ``-m ea`` for end and abortion messages.

   * - ``-a XX``
     - ``--begin=XX``
     - Specify the earliest time for a job to start.
       SGE format:  ``[YYMMDDhhmm]``.
       Slurm format: ``YYYY-MM-DD[HH:MM[:SS]]``.

   * - ``-wd XX``
     - ``--workdir=XX``
     - Execute the job from the specified directory.

   * - ``-l excl=true``
     - ``--exclusive``
     - Request exclusive access to all nodes used by the job so no other jobs can
       run on them.  This can be useful for benchmarking purposes where you want
       to ensure that you have exclusive use of e.g. memory/IO buses.  Note that
       you still need to request CPU cores and memory to avoid being limited to
       just the default per job (one core and a set amount of RAM).  Also note
       that the use of this option will likely result in longer queuing times.

   * - ``-P XX`` (and possibly ``-q YY``)
     - ``-p XX`` and ``-q YY``
     - Request a non-standard session, typically on private nodes
       and/or nodes that have partiticular hardware by,
       on SGE, requesting a Project and possibly a Queue,
       or, on Slurm, requesting a Partition and optionally a QoS.
       You should have been told exactly what to request
       when you were granted access to particular hardware.

   * - ``-l hostname=XX``
     - ``--nodelist=XX``
     - Target a node by name.
       Not recommended for normal use.

To learn more about all available options run one of the following whilst on a relevant cluster:

* ``man qsub``
* ``man sbatch``

Converting SGE job scripts to Slurm scripts
-------------------------------------------

Stanford Research Computing Center provide a `guide to this <https://srcc.stanford.edu/sge-slurm-conversion>`_.
