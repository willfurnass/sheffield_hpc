.. _dgx1_com_groupnodes_sharc:

Nvidia DGX-1 [COM]
==================

   The Nvidia DGX-1 is the world's first Deep Learning supercomputer.
   It is equipped with 8 Tesla P100 GPUs
   connected together with NVLink technology to provide super fast inter-GPU communication.
   Capable of performing 170 TeraFLOPs of computation,
   it can provide up to 75 times speed up on training deep neural networks 
   compared to the latest Intel Xeon CPU.

DGX-1 Specifications
--------------------

* 8x Tesla P100 GPUs (16GB RAM each)
* Dual 20-core Intel Xeon E5-2698 v4 2.2 GHz
* 512 GB System RAM

.. note::

   One GPU is faulty; only 7 GPUs are currently usable.

Requesting Access
-----------------

Same as per the :ref:`COM big memory nodes <big_mem_com_groupnodes_sharc>`.

Starting an interactive session
-------------------------------

Once you have been granted access to the RSE nodes in ShARC, 
to request an interactive session (interactive job) on the DGX-1 node with a single GPU, 
type:

.. code-block:: sh

   qrshx -l gpu=1 -P rse -q rse-interactive.q

* ``-l gpu=1`` denotes the number of GPUs that will be used in the job (maximum of 7), 
  This is required otherwise the job will be placed on a one of the RSE team's CPU-only nodes.
* ``-P rse -q rse.q`` denotes that the job will be submitted under the ``rse`` Project and 
  you want it to only run in the ``rse-interactive.q`` job queue,
  which is for interactive jobs that can run for up to 8 hours.

.. note::

   You are limited to at most 1 GPU over all your jobs on the DGX-1 that do not run in the ``rse.q`` (batch-job-only) job queue.
   This is to encourage users to prefer batch jobs over interactive sessions on the DGX-1, 
   as batch jobs make more efficient use of resources.

Submitting batch jobs
---------------------

Batch jobs can be submitted to the DGX-1 by 
adding lines containing ``-l gpu``, ``-P rse`` and ``-q rse.q`` 
to your job submission script (note the different argument to ``-q``). 

For example, create a job script named ``my_job_script.sh`` with the contents:

.. code-block:: sh

   #!/bin/bash
   #$ -l gpu=1 
   #$ -P rse 
   #$ -q rse.q

   echo "Hello world"

You can add additional lines beneath ``#$ -q rse.q`` to request additional resources 
e.g. ``-l rmem=10G``  to request 10GB RAM per CPU core rather than the default.

Run your script with the ``qsub`` command:

.. code-block:: sh

   qsub my_job_script.sh

You can use ``qstat`` command to check the status of your current job. 
An output file is created in your home directory that captures your script's outputs.

Note that the maximum run-time for jobs submitted to the (batch job only) ``rse.q`` is four days, 
as is standard for batch jobs on ShARC.

See :ref:`submit-queue` for more information on job submission and the Sun Grid Engine scheduler.

Deep Learning on the DGX-1
--------------------------

Many popular Deep Learning packages are available to use on the DGX-1 and the ShARC cluster. 
Please see :ref:`DeepLearning_sharc` for more information.
