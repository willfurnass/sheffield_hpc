ANSYS
=====

.. sidebar:: ANSYS
   
   :Versions: 19.4, 20.1, 20.2 
   :Dependencies: UDFs / User subroutines need the GCC/7.3.0-2.30 compiler.
   :URL: http://www.ansys.com 

The ANSYS suite of programs can be used to numerically simulate a large variety of structural and fluid dynamics problems found in many engineering, physics, medical, aeronotics and automative industry applications.

Interactive usage
-----------------

After connecting to Bessemer (see :ref:`ssh`),  start an `interactive graphical session <https://docs.hpc.shef.ac.uk/en/latest/hpc/scheduler/submit.html#interactive-sessions>`_.

ANSYS version 20.2 can be activated using the module file::

    module load ANSYS/20.2/binary

and the workbench is launched using::

    runwb2

Fluent, CFX, ICEM, Mechanical APDL/model (and many more) can all be accessed from the workbench. Outside the workbench the corresponding GUIs can be launched using ``fluent``, ``cfx5pre``, ``icemcfd`` and ``launcher``.

ANSYS example models
--------------------

ANSYS contains a large number of example models which can be used to become familiar with the software.
The models can be found in::

    /usr/local/packages/live/noeb/ANSYS/20.2/binary/v202/ansys/data/
	

Batch jobs
----------

The following is an example batch submission script, ``cfd_job.sh``, to run the executable ``fluent`` with input journal file ``subjou.jou``, and carry out a 2D double precision CFD simulation. The script requests 4 cores using the OpenMP parallel environment with a runtime of 60 mins and 2 GB of real memory per core. ::

    #!/bin/bash
    #SBATCH --nodes=1
    #SBATCH --ntasks-per-node=4
    #SBATCH --mem=8000
    #SBATCH --job-name=name_fluent_smp_4
    #SBATCH --output=output_fluent_smp_4
    #SBATCH --time=01:00:00
    #SBATCH --mail-user=joe.bloggs@sheffield.ac.uk
    #SBATCH --mail-type=ALL
    module load ANSYS/20.2
    fluent 2ddp -i subjou.jou -g -t4
	
The job is submitted to the queue by typing::

    sbatch cfd_job.sh

The following is an example batch submission script, ``mech_job.sh``, to run the mechanical executable ``mapdl`` with input file ``CrankSlot_Flexible.inp``, and carry out a mechanical simulation. The script requests 2 cores using the OpenMP parallel environment with a runtime of 60 mins and 2 GB of real memory per core. ::

    #!/bin/bash
    #SBATCH --nodes=1
    #SBATCH --ntasks-per-node=2
    #SBATCH --mem=4000
    #SBATCH --job-name=ansys_mech-test
    #SBATCH --output=output_ansys_mech_test
    #SBATCH --time=01:00:00
    #SBATCH --mail-user=joe.bloggs@sheffield.ac.uk
    #SBATCH --mail-type=ALL
    module load ANSYS/20.2
    ANSYS_OPTIONS="-smp -dir $(pwd) -b -np $SLURM_NTASKS -j solution -i"
    mapdl $ANSYS_OPTIONS CrankSlot_Flexible.inp

The job is submitted to the queue by typing::

    sbatch mech_job.sh

Installation note for Administrators:
-------------------------------------

mapdl will not run without modifying the file::

    /usr/local/packages/live/noeb/ANSYS/20.2/binary/v202/ansys/bin/anssh.ini

The following instruction should be inserted at line 2433 in ``anssh.ini``::

    setenv KMP_AFFINITY compact

	
Please note ANSYS 20.1 and 20.2 have been installed manually with the GUI in the following directories and permissions corrected as follows::
	
    /usr/local/packages/live/noeb/ANSYS/20.1/binary/
    /usr/local/packages/live/noeb/ANSYS/20.2/binary/
	
    chmod 775 -R /usr/local/packages/live/noeb/ANSYS/20.1/binary/
    chmod 775 -R /usr/local/packages/live/noeb/ANSYS/20.2/binary/
	
Please follow the same install directory structure.

In addition the following software packages are not included with the installations::


    "ANSYS Chemkin"
    "ANSYS Geometry Interfaces".
	

