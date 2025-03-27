The necessary scripts for the upcoming exercises are located in our
`hpc-examples <https://github.com/rcgsheffield/hpc-examples>`_ repository.
This repository is accessible on our Stanage HPC cluster. To utilise it, load the module:

.. code-block:: bash

   module load hpc-examples

After loading, you can access the tutorial scripts via the ``$HPC_TUTS`` environment variable.

For example, you can then copy ``slurm/pi.py`` to your current working directory:

.. code-block:: bash

   cp $HPC_TUTS/slurm/pi.py . 
