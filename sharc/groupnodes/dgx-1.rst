.. _dgx1_com_groupnodes_sharc:

Nvidia DGX-1 [COM]
==================

The Nvidia DGX-1 is the world's first Deep Learning supercomputer. It is equipped with 8 Tesla P100 GPUs connected together with NVLink technology to provide super fast inter-GPU communication. Capable of performing 170 Teraflops of computation, it can provide up to 75 times speed up on training deep neural networks compared to the latest Intel Xeon CPU.

DGX-1 Specifications
--------------------

* 8x Tesla P100 GPUs (16GB RAM each)
* Dual 20-core Intel Xeon E5-2698 v4 2.2 GHz
* 512 GB System RAM

Requesting Access
-----------------

The node managed by the `RSE group <https://rse.shef.ac.uk>`_ and is available by request to all research groups belonging to the Computer Science Department.

To use the DGX-1 you must join the computer science RSE group and submit jobs to the RSE queue. To join this group please email `Twin Karmakharm <t.karmakharm@sheffield.ac.uk>`_  or `RSE enquiries <rse@shef.ac.uk>`_.

Starting an interactive session
-------------------------------

Once you have obtained permission to use the node, to request an interactive DGX-1 node, type: ::

	qrshx -l gpu=1 -P rse -q rse-interactive.q

You will then be placed in a special RSE queue that has the DGX-1. The parameter ``-l gpu=1`` determines the number of GPUs that will be used in the job (maximum of 8), ``-P rse -q rse.q`` then denotes that you run the job under the ``rse`` project with the ``rse-interactive.q`` queue (which is for interactive jobs that can run for up to 8 hours).  Note that the parameter ``-l gpu=`` is required otherwise you will be placed on a one of the RSE team's CPU-only nodes.

Submitting batch jobs
---------------------

Batch Jobs can be submitted to the DGX-1 by adding the ``-l gpu=1 -P rse -q rse.q`` parameters (note the different argument to ``-q``). 
For example, create a job script named ``my_job_script.sh`` with the contents: ::

	#!/bin/bash
	#$ -l gpu=1 
        #$ -P rse 
        #$ -q rse.q

	echo "Hello world"

You can add additional lines beneath ``#$ -q rse.q`` to request addtional resources e.g. ``-l rmem=10G``  to request 10GB RAM per CPU core rather than the default.

Run your script with the ``qsub`` command ::

	qsub my_job_script.sh

You can use ``qstat`` command to check the status of your current job. An output file is created in your home directory that captures your script's outputs.

Note that the maximum run-time for jobs submitted to the (batch job only) ``rse.q`` is four days, as is standard for batch jobs on ShARC.

See :ref:`submit-queue` for more information on job submission and the Sun Grid Engine scheduler.

Deep Learning on the DGX-1
--------------------------

Many popular Deep Learning packages are available to use on the DGX-1 and the ShARC cluster. Please see :ref:`DeepLearning_sharc` for more information.
