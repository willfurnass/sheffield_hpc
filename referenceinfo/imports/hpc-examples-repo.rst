.. note::

     The necessary scripts for the upcoming exercises are located in our
     `hpc-examples <https://github.com/rcgsheffield/hpc-examples>`_ repository.
     This repository is accessible on our Stanage HPC cluster. To utilise it, load the module:
     
     .. code-block:: bash
     
        module load hpc-examples
     
     After loading, you can access the examples scripts via the ``$HPC_EXAMPLES`` environment variable.
     
     For example, you can then run ``slurm/pi.py`` in the following way:
     
     .. code-block:: bash
     
        python $HPC_EXAMPLES/slurm/pi.py  
