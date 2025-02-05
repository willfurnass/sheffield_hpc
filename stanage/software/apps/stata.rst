.. _stata_stanage:

.. |softwarename| replace:: STaTa
.. |currentver| replace:: 18

STaTa
======

.. warning::
    The stata site licence only allows 4 consecutive cores per job. Meaning batch jobs cant use more than 4 cores. Please *DO NOT* request more than 4 cores as the additional cores will not be used, lowering your efficiency. 


.. include:: /referenceinfo/imports/stanage/packages/stata-sdbr-el7-icelake-znver-stanage.rst

.. include:: /referenceinfo/imports/stanage/packages/stata-dscr-el7-icelake-znver-stanage.rst


--------

Interactive usage
-----------------

.. include:: /referenceinfo/imports/scheduler/SLURM/common_commands/srun_start_interactive_session_import_stanage.rst

STaTa can be loaded with the command:

.. include:: /referenceinfo/imports/stanage/packages/stata-ml-el7-icelake-znver-stanage.rst

After this any of the STAR commands can be run from the terminal prompt. The available 
commands can be obtained using:

.. code-block:: bash

	stata -h

--------

Batch usage
-----------

The following is an example batch submission script, ``my_job.sh``, to run the executable ``stata``.
The script requests 4 cores with a runtime of 5 minutes and 1 GB of memory.

Stata script ``hello_world.do`` : 

.. code-block:: bash

        disp "Hello world" 

Batch script ``my_job.sh`` :

.. code-block:: bash

        #!/bin/bash
        #SBATCH --job-name=stata_test
        #SBATCH --cpus-per-task=4
        #SBATCH --mem=1000
        #SBATCH --output=output_stata_4.%j.out
        #SBATCH --time=00:05:00
        #SBATCH --mail-user=a.person@sheffield.ac.uk
        #SBATCH --mail-type=ALL

        module load stata/18.0
        stata -b hello_world.do      

The job is submitted to the queue by typing:

.. code-block:: console

   $ sbatch my_job.sh



--------

Installation notes
------------------

Installation method
^^^^^^^^^^^^^^^^^^^

This section is primarily for administrators of the system. |softwarename| has been installed using the default Easybuild config files.

Build logs and test reports can be found in ``$EBROOTSTATA/easybuild`` with a given module loaded.

The stata license was added manually by following the instructions `here <https://www.stata.com/install-guide/change-information/>`_.

--------

Testing
^^^^^^^
Testing was performed using the above example.

.. include:: /referenceinfo/imports/stanage/packages/stata-dpnd-el7-icelake-znver-stanage.rst

