.. _big_mem_com_groupnodes_sharc:

Big memory nodes [COM]
======================

The Department of Computer Science have purchased 8 nodes in ShARC that each have 
much more than the standard amount of RAM per node. 

Specifications
--------------

* 3x nodes with
   * 2x Intel Xeon E5-2630 v3 processors (2.40 GHz, 8 cores per socket i.e. 16 total)
   * 768 GB RAM (48 GB / CPU core)
* 5x nodes with
   * 2x Intel Xeon Gold 6138 processors (2.00GHz, 20 cores per socket i.e. 40 total)
   * 768 GB RAM (19.2 GB / CPU core)
   * 1TB SSD

Requesting Access
-----------------

The nodes are managed by the `RSE group <http://rse.shef.ac.uk>`_ and are available by request to all research groups belonging to the Computer Science Department.

To use the nodes you must 

#. Be made a member of the ``rse`` Grid Engine (scheduler) *Access Control List* (ACL i.e. user group);
#. Submit jobs using the ``rse`` Grid Engine *Project*;
#. Start interactive jobs in ``rse-interactive.q`` Grid Engine *Cluster Queue*;
#. Start batch jobs in the ``rse.q`` Grid Engine *Cluster Queue*;

To join this ACL please email `RSE enquiries <rse@shef.ac.uk>`_.

Running an interactive session
------------------------------

Once you have obtained permission to use the nodes you can request an interactive session on one of the nodes using:

.. code-block:: bash

	qrshx -P rse -q rse-interactive.q

Here ``-P rse`` specifies that you want to use the ``rse`` project for your session, 
which gives you access to these big memory nodes and 
ensures that your interactive session can run in the ``rse-interactive.q`` job queue 
(as can be seen if you subsequently run ``qstat -u $USER`` from within your session).

The ``rse-interactive.q`` job queue has a maximum job runtime (``h_rt``) of four hours, 
as is standard for interactive jobs on SHARC.

Submitting batch jobs
---------------------

Jobs can be submitted to the nodes by adding the ``-P rse`` and ``-q rse.q`` parameters. 
For example, create a job script named ``my_job_script.sh`` with the contents:

.. code-block:: bash

	#!/bin/bash
	#$ -P rse 
	#$ -q rse.q

	echo "Hello world"

You can of course add more options to the script such as a request for additional RAM (e.g. ``$# -l rmem=10G``).

Run your script with the ``qsub`` command:

.. code-block:: bash

	qsub my_job_script.sh

You can use the ``qstat`` command to check the status of your current job. 
An output file is created in your home directory that captures your script's outputs.

See :ref:`submit-queue` for more information on job submission and the Sun Grid Engine scheduler.
