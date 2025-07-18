.. _srun:

srun
====

``srun`` is the SLURM command used to launch tasks within a job allocation. It can be used in both batch scripts and on the command line to:

* Start one or more job 'steps' in a batch job (a step can request and use some or all of the resources allocated to the job; steps can run sequentially or in parallel) 
    * useful for analysing performance across multiple stages of a workflow
* Launch :ref:`MPI <parallel_MPI>` or multi-task jobs
    * (e.g. with ``--ntasks`` > 1)
    * where multiple processes need to be started across CPUs or nodes
* Run an interactive job 
    * for example, ``srun --pty bash -i``

Documentation
-------------

Run ``man srun`` on the system for full details.

Usage
-----

To launch an interactive shell:

.. code-block:: slurm

    srun --pty bash -i

To run a non-interactive job directly:

.. code-block:: slurm

    srun --export=ALL --time=00:10:00 --mem=1G python3 my_script.py 100000000

In a batch script:

.. code-block:: slurm

    #SBATCH --ntasks=4
    srun --export=ALL ./my_parallel_app

For simple, single-task jobs, ``srun`` is optional in batch scripts — your command can be run directly.

See the :ref:`interactive and batch submission guide <submit_interactive_stanage>` for more examples.
