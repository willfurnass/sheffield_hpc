.. _dcs_gpu_nodes_bessemer:

GPU nodes for all Computer Science researchers
==============================================

GPU nodes purchased for :ref:`Bessemer <bessemer>` by the `Department of Computer Science <https://www.sheffield.ac.uk/dcs>`__ (DCS)
for use by DCS research staff, their collaborators and their research students.

.. _dcs_gpu_nodes_hw:

Hardware specifications
-----------------------

Eight nodes (``bessemer-node030`` to ``bessemer-node037``) each have:

.. list-table::
   :header-rows: 0

   * - Processors
     - 2x Intel Xeon Gold 6138 (2.00GHz; 40 cores per CPU)
   * - RAM
     - 192GB (DDR4 @ 2666 MHz)
   * - NUMA nodes
     - 2x
   * - GPUS
     - 4x NVIDIA Tesla V100 SXM2 (32GB RAM each; NVLINK interconnects between GPUs)
   * - Networking
     - 25 Gbps Ethernet
   * - Local storage
     - 140 GB of temporary storage under ``/scratch`` (2x SSD RAID1)

Requesting access
-----------------

Contact ``rse@sheffield.ac.uk``.

.. todo::
   Update this section

Using the nodes
---------------

There are several ways to access these nodes.
The type of access granted for a job depends on which Slurm *Account* and *Partition* are requested at job submission time.

1. DCS test/debugging access
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

E.g. for short test batch jobs or for interactive debugging.

To submit a job via this route, you need to :ref:`specify a *Partition* and *Account* <slurm_access_priv_nodes>` when submitting a batch job or starting an interactive session:

* Partition: ``dcs-gpu-test``
* Account: ``dcs-gpu``

Resource limits per job:

* Exactly 1 or 2 GPUs must be requested
* Default run-time: 30 minutes
* Maximum run-time: 30 minutes
* :ref:`Number of CPU cores, amount of RAM and number of GPUs in a single node <dcs_gpu_nodes_hw>`
  i.e. multi-node jobs are not permitted.

Each user can run a maximum of two of these jobs concurrently.

2. DCS access for larger jobs
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

If you want to run a longer job that uses up to *all* the resources available in *one* of these nodes
then you can :ref:`specify a different Partition <slurm_access_priv_nodes>` when submitting a batch job or starting an interactive session:

* Partition: ``dcs-gpu``
* Account: ``dcs-gpu``

Please *only run batch jobs this way*: long-running interactive sessions that are associated with large resource requests are often an inefficient way of using cluster resources.

Resource limits per job:

* At least one GPU must be requested
* Default run-time: 8 hours
* Maximum run-time: 7 days
* :ref:`Number of CPU cores, amount of RAM and number of GPUs in a single node <dcs_gpu_nodes_hw>`
  i.e. multi-node jobs are not permitted.

.. todo::
   If using cluster-wide values for default and max run time then link to central info re that rather than duplicating here.

.. todo::

   Leave commented until implemented and tested

   3. General pre-emptable access

   Users other than Computer Science researchers and their collaborators can
   make use of idle time on these nodes and other nodes
   for running GPU jobs *or* CPU-only jobs
   by submitting batch jobs and starting interactive sessions :ref:`specifying a particular partition <slurm_access_priv_nodes>`:

   * Partition: ``preempt``

   These jobs can be pre-empted by jobs submitted to the ``dcs-gpu`` and ``dcs-gpu-test`` partitions;
   if this happens
   the pre-empted jobs will be stopped mid-execution and re-queued.
