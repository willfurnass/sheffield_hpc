
Submitting Non-Interactive Jobs
-------------------------------

Additional options for job submission
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Specify nodes and tasks for MPI jobs:

.. code-block:: sh

    #SBATCH --nodes=1
    #SBATCH --ntasks-per-node=16

Memory allocation:

.. code-block:: sh

    #SBATCH --mem=16000

Specify the output file name:

.. code-block:: sh

    #SBATCH --output=output.%j.test.out

Key Slurm Scheduler Commands
----------------------------

Display the job queue. Jobs typically pass through several states in the course of their execution. The typical states are PENDING, RUNNING, SUSPENDED, COMPLETING, and COMPLETED.

.. code-block:: sh

    squeue

Shows job details:

.. code-block:: sh

    sacct -v

Details the HPC nodes:

.. code-block:: sh

    sinfo
