Abaqus
======

.. sidebar:: Abaqus
   
   :Versions: 6.14.2 (see Addendum section), 2018 
   :Dependencies: User subroutines need the Intel FORTRAN compiler 2019.
   :URL: http://www.3ds.com/products-services/simulia/products/abaqus/ 
   :Documentation: https://www.3ds.com/products-services/simulia/products/abaqus/ (note: register for an account to access)

Abaqus is a software suite for Finite Element Analysis (FEA) developed by Dassault Systèmes.

Interactive usage
-----------------

After connecting to Bessemer (see :ref:`ssh`),  start an `interactive graphical session <https://docs.hpc.shef.ac.uk/en/latest/hpc/scheduler/submit.html#interactive-sessions>`_.

Abaqus version 2018 can be activated using the module file::

    module load ABAQUS/2018

and launched using::

    abaqus cae


**Note:** there is an Abaqus/SLURM mpi issue which prevents Abaqus from running jobs correctly. To rectify this problem the following command should be used prior to lauching the GUI::

    unset SLURM_GTIDS

------------

Abaqus example problems
-----------------------

Abaqus contains a large number of example problems which can be used to become familiar with Abaqus on the system.
These example problems are described in the Abaqus documentation and can be obtained using the Abaqus ``fetch`` command.
For example, after loading the Abaqus module enter the following at the command line to extract the input file for test problem s4d::

    abaqus fetch job=s4d
	
This will extract the input file ``s4d.inp`` to run the computation defined by the commands and batch submission script below.

------------

Batch jobs
----------

The following is an example batch submission script, ``my_job.sh``, to run the executable ``abaqus`` with input file is ``s4d.inp``. The script requests 4 cores using the OpenMP parallel environment ``smp`` with a runtime of 30 mins and 2 GB of real memory per core. ::

    #!/bin/bash
    #SBATCH --comment=abaqus_smp_test
    #SBATCH --nodes=1
    #SBATCH --ntasks-per-node=4
    #SBATCH --mem=8000
    #SBATCH --output=output_abaqus_smp_4
    #SBATCH --time=00:30:00
    #SBATCH --mail-user=joe.bloggs@sheffield.ac.uk
    #SBATCH --mail-type=ALL
    module load ABAQUS/2018
    unset SLURM_GTIDS
    abaqus job=my_job input=s4d.inp mp_mode=threads cpus=$SLURM_NTASKS scratch=$TMPDIR memory="8gb" interactive
	
The job is submitted to the queue by typing::

    sbatch my_job.sh

**User subroutines:** The ``umatmst3`` model has a user defined subroutine ``umatmst3.f``. The model files are obtained using ::

    abaqus fetch job=umatst3*

The script below is an example of a batch submission script for a single core job with a runtime of 30 mins, 8 GB of real memory and with user subroutine ``umatmst3.f`` and input file ``umatmst3.inp``. ::

    #!/bin/bash
    #SBATCH --comment=abaqus_subroutine_test
    #SBATCH --nodes=1
    #SBATCH --ntasks-per-node=1
    #SBATCH --mem=8000
    #SBATCH --output=output_abaqus_subroutine
    #SBATCH --time=00:30:00
    #SBATCH --mail-user=joe.bloggs@sheffield.ac.uk
    #SBATCH --mail-type=ALL
    module load ABAQUS/2018
    module load ifort/2019.1.144-GCC-8.2.0-2.31.1
    unset SLURM_GTIDS
    abaqus job=my_job input=umatmst3.inp user=umatmst3.f scratch=$TMPDIR memory="8gb" interactive

Note that the module ``ifort/2019.1.144-GCC-8.2.0-2.31.1``, required for compiling the user subroutines, is not automatically loaded when the module for Abaqus is loaded.

------------

Licensed options
----------------

All available Abaqus licenses can be viewed using ``abaqus licensing r`` e.g. ::

   $ module load ABAQUS/2018
   $ abaqus licensing r

   Feature                         Version     #licenses    Expires      Vendor
   _______                         _________   _________    __________   ______
   abaqus_extended                 61.9         19          31-dec-2018  ABAQUSLM
   abaqus                          61.9         250         31-dec-2018  ABAQUSLM
   ams                             61.9         1           31-dec-2018  ABAQUSLM
   aqua                            61.9         250         31-dec-2018  ABAQUSLM
   cosim_acusolve                  61.9         1           31-dec-2018  ABAQUSLM
   cosim_direct                    61.9         1           31-dec-2018  ABAQUSLM
   cse                             61.9         1           31-dec-2018  ABAQUSLM
   design                          61.9         250         31-dec-2018  ABAQUSLM
   euler_lagrange                  61.9         1           31-dec-2018  ABAQUSLM
   gpgpu                           61.9         1           31-dec-2018  ABAQUSLM
   multiphysics                    61.9         1           31-dec-2018  ABAQUSLM
   parallel                        61.9         16384       31-dec-2018  ABAQUSLM
   sw_assoc_import                 61.9         1           31-dec-2018  ABAQUSLM
   catiav5_assoc_import            61.9         1           31-dec-2018  ABAQUSLM
   catiav5_import                  61.9         1           31-dec-2018  ABAQUSLM
   catiav6_assoc_import            61.9         1           31-dec-2018  ABAQUSLM
   tomee                           61.9         1           31-dec-2018  ABAQUSLM
   pydriver                        61.9         1           31-dec-2018  ABAQUSLM
   cae                             61.9         19          31-dec-2018  ABAQUSLM
   rtgateway                       61.9         19          31-dec-2018  ABAQUSLM
   gateway                         61.9         19          31-dec-2018  ABAQUSLM
   safe_ex_gui                     61.9         19          31-dec-2018  ABAQUSLM
   cfd                             61.9         250         31-dec-2018  ABAQUSLM
   explicit                        61.9         250         31-dec-2018  ABAQUSLM
   foundation                      61.9         250         31-dec-2018  ABAQUSLM
   simflow                         61.9         250         31-dec-2018  ABAQUSLM
   standard                        61.9         250         31-dec-2018  ABAQUSLM
   cse_token                       61.9         250         31-dec-2018  ABAQUSLM
   safe_ex_engine                  61.9         250         31-dec-2018  ABAQUSLM
   tosca_topo                      61.9         250         31-dec-2018  ABAQUSLM
   tosca_shape                     61.9         250         31-dec-2018  ABAQUSLM
   tosca_bead                      61.9         250         31-dec-2018  ABAQUSLM
   tosca_sizing                    61.9         250         31-dec-2018  ABAQUSLM
   tosca_int_abaqus                61.9         250         31-dec-2018  ABAQUSLM
   tosca_int_ansys                 61.9         250         31-dec-2018  ABAQUSLM
   tosca_int_nastran               61.9         250         31-dec-2018  ABAQUSLM
   tosca_adv_nonlinear             61.9         250         31-dec-2018  ABAQUSLM
   tosca_adv_durability            61.9         250         31-dec-2018  ABAQUSLM
   tosca_adv_morph                 61.9         250         31-dec-2018  ABAQUSLM
   tosca_smooth                    61.9         250         31-dec-2018  ABAQUSLM
   tosca_report                    61.9         250         31-dec-2018  ABAQUSLM
   tfluid_topo                     61.9         250         31-dec-2018  ABAQUSLM
   tfluid_smooth                   61.9         250         31-dec-2018  ABAQUSLM
   tfluid_parallel                 61.9         250         31-dec-2018  ABAQUSLM
   tfluid_int_ccmp                 61.9         250         31-dec-2018  ABAQUSLM
   tfluid_int_fluent               61.9         250         31-dec-2018  ABAQUSLM

Run ``abaqus licensing`` for usage info for the Abaqus licensing sub-command. Run ``abaqus licensing ru`` to see current licence usage.

------------

Checkpointing your work
-----------------------

Abaqus has a built-in checkpoint and restart feature.

Add the following to the input file (refer to official Abaqus documentation for detail): ::

   *RESTART, WRITE, OVERLAY, FREQUENCY=10

**OVERLAY** saves only one state, i.e. overwrites the restart file every time new restart information is written
    
**FREQUENCY=N** writes restart information every N timesteps

And, to restart the job, create a new input file newJobName with only a single line:  ::

   *RESTART, READ

Then run Abaqus specifying both the new and old job names:  ::

   abaqus jobname=newJobName oldjob=oldJobName

------------

Addendum: Abaqus 6.14.2 (non-EasyBuild install):
------------------------------------------------

Abaqus 6.14.2 was installed using the standard Abaqus installer due to issues using EasyBuild.

It can be activated using the following module commands::

    module use /usr/local/modulefiles/live/apps
    module load ABAQUS/6.14.2/binary

and launched using::

    abaqus cae

The following is an example batch submission script, ``my_job.sh``, to run the executable ``abaqus`` with input file is ``s4d.inp``. The script requests 4 cores using the OpenMP parallel environment ``smp`` with a runtime of 30 mins and 2 GB of real memory per core. ::

    #!/bin/bash
    #SBATCH --comment=abaqus_smp_test
    #SBATCH --nodes=1
    #SBATCH --ntasks-per-node=4
    #SBATCH --mem=8000
    #SBATCH --output=output_abaqus_smp_4
    #SBATCH --time=00:30:00
    #SBATCH --mail-user=joe.bloggs@sheffield.ac.uk
    #SBATCH --mail-type=ALL
    module use /usr/local/modulefiles/live/apps
    module load ABAQUS/6.14.2/binary
    unset SLURM_GTIDS
    abaqus job=my_job input=s4d.inp mp_mode=threads cpus=$SLURM_NTASKS scratch=$TMPDIR memory="8gb" interactive

