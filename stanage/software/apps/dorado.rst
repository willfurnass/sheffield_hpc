.. |softwarename| replace:: Dorado
.. |currentver| replace:: 1.0.1

.. _dorado_stanage:

Dorado
======

.. sidebar::  |softwarename|

   :Versions: |currentver|
   :Dependencies: Apptainer
   :URL: https://dorado-docs.readthedocs.io/en/latest/
   :Documentation: https://dorado-docs.readthedocs.io/en/latest/

Dorado is a high-performance, easy-to-use, open source analysis engine for Oxford Nanopore reads.

.. warning:: 

   Dorado requires a dependencies not available on CentOS 7. To make it run on stanage the team containerized the installation using Apptainer. As a result, **Dorado is only recommended for use on single-node jobs**. The container integration has been abstracted via the module system, allowing you to use Dorado as if it were a standard module.   
   
   **Note: This is an experimental build. Please report any issues to the HPC support team.**

Interactive Usage
-----------------
.. include:: /referenceinfo/imports/scheduler/SLURM/common_commands/srun_start_interactive_session_import_stanage.rst


.. code-block:: bash

    #This will load the Dorado module and make the command available in your path.
    module load dorado/1.0.1
    
    # Run Dorado basecaller
    dorado basecaller hac pod5s/ > calls.bam

Batch Submission
----------------


.. code-block:: bash

    #!/bin/bash
    #SBATCH --nodes=1
    #SBATCH --mem=5G
    #SBATCH --time=00:30:00
    #SBATCH --cpus-per-task=4
    #Email notifications to me@somedomain.com
    #SBATCH --mail-user=me@somedomain.com
    #Email notifications if the job fails
    #SBATCH --mail-type=FAIL
    #SBATCH --job-name=dorado_test
    
    module load dorado/1.0.1
    
    dorado basecaller hac pod5s/ > calls.bam

Installation Notes
------------------
Dorado requires a newer version of GLIBC not available on the CentOS 7, so it is installed using Apptainer. The container definition file and easyconfig can be found `here <https://github.com/rcgsheffield/hpc-apptainer-images/tree/main/software/dorado/1.0.1>`_ 


Testing
-------


