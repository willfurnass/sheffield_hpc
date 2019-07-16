Abaqus
======

.. sidebar:: Abaqus
   
   :Versions: 2018
   :URL: http://www.3ds.com/products-services/simulia/products/abaqus/ 
   :Documentation: https://www.3ds.com/support/documentation/users-guides/
   :Local URL: https://www.sheffield.ac.uk/cics/research/software/abaqus


Abaqus is a software suite for Finite Element Analysis (FEA) developed by Dassault Systèmes.


Usage
-----

Abaqus versions 2018 can be activated using the module files::

   module use /usr/local/modulefiles/staging/eb/cae/
   module load ABAQUS/2018
	
Type ``abaqus cae`` to launch the Abaqus GUI from an interactive session with X Window support. Please see usage note below for graphics support options.
Type ``abaqus`` for the command line interface. Typing ``abaqus -help`` will display a list of usage options.

Batch jobs
----------

The easiest way of running a batch job for a particular version of Abaqus is::
    
    module use /usr/local/modulefiles/staging/eb/cae/
    module load ABAQUS/2018
    abaqus input=static_square.inp

Users are encouraged to write their own batch submission scripts. The following is an example batch submission script, ``my_job.sh``, to run the executable ``abaqus`` and which is submitted to the queue by typing ``sbatch my_job.sh``::

    #!/bin/bash
    #SBATCH --comment=abaqus_test
    #SBATCH --nodes=1
    #SBATCH --ntasks-per-node=16
    #SBATCH --mem=64000
    #SBATCH --output=output.%j.abaqus.test
    #SBATCH --time=00:30:00
    #SBATCH --mail-user=$USER@sheffield.ac.uk

    module use /usr/local/modulefiles/staging/eb/cae/
    module load ABAQUS/2018

    abaqus input=static_square.inp \
        job=OUT_static_square \
        cpus=$SLURM_NTASKS \
        mp_mode=THREADS
	
The above script requests 16 cores with a runtime of 30 mins and 64 GB of real memory per core. The Abaqus input file is ``static_square.inp``.
