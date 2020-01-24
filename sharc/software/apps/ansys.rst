ANSYS
=====

.. sidebar:: ANSYS
   
   :Versions: 15.0, 16.1, 17.2, 18.0, 18.2, 19.0, 19.1, 19.2, 19.3 & 19.4
   :Dependencies: No prerequsite modules loaded. However, If using the User Defined Functions (UDF) will also need the following: For ANSYS Mechanical, Workbench, CFX and AutoDYN: Intel 14.0 or above; Compiler For Fluent: GCC 4.6.1 or above
   :URL: http://www.ansys.com 
   :Local URL: https://www.shef.ac.uk/cics/research/software/fluent


The Ansys suite of programs can be used to numerically simulate a large variety of structural and fluid dynamics problems found in many engineering, physics, medical, aeronotics and automative industry applications.


Usage
-----

Ansys can be activated using the module files::

    module load apps/ansys/15.0
    module load apps/ansys/16.1
    module load apps/ansys/17.2
    module load apps/ansys/18.0/binary
    module load apps/ansys/18.2/binary
    module load apps/ansys/19.0/binary
    module load apps/ansys/19.1/binary
    module load apps/ansys/19.2/binary
    module load apps/ansys/19.3/binary
    module load apps/ansys/19.4/binary
	

The Ansys Workbench GUI executable is ``ansyswb``. ``ansyswb`` can be launched during an interactive session with X Window support (e.g. an interactive ``qrshx`` session).
The Ansys Mechanical executable is ``ansys-mechanical`` and ``fluent`` is the executable for Fluent.
 
NOTE: for Ansys versions >= 18.0 using the command ``fluent`` results in an unresponsive fluent launcher. To launch fluent and bypass the launcher use ``fluent dim`` where dim = 2d, 2ddp, 3d or 3ddp.

ANSYS example models
--------------------

ANSYS contains a large number of example models which can be used to become familiar with the software.
The models can be found in::

    /usr/local/packages/live/eb/ANSYS/19.4/v194/ansys/data
	

Batch jobs
----------
	
``Fluent CFD``: the following is an example batch submission script, ``cfd_job.sh``, to run the executable ``fluent`` with input journal file ``test.jou``, and carry out a 2D double precision CFD simulation. The script requests 8 cores using the MPI parallel environment ``mpi-rsh`` with a runtime of 30 mins and 2 GB of real memory per core. The Fluent input journal file is ``test.jou``. **Note:** Please use the ``mpi-rsh`` parallel environment to run MPI parallel jobs for Ansys/Fluent. ::

    #!/bin/bash
    #$ -cwd
    #$ -M joe.bloggs@sheffield.ac.uk
    #$ -m abe
    #$ -l h_rt=00:30:00
    #$ -l rmem=2G
    #$ -pe mpi-rsh 8

    module load apps/ansys/19.4

    fluent 2ddp -i test.jou -g -t$NSLOTS -sge -mpi=intel -rsh -sgepe mpi-rsh

The job is submitted to the queue by typing::

    qsub cfd_job.sh

``Mapdl mechanical``: the following is an example batch submission script, ``mech_job.sh``, to run the mechanical executable ``mapdl`` with input file ``CrankSlot_Flexible.inp``, and carry out a mechanical simulation. The script requests 4 cores using the OpenMP (``single node shared memory``) parallel environment with a runtime of 10 mins and 2 GB of real memory per core. ::

    #!/bin/bash
    #$ -cwd
    #$ -N mech_test
    #$ -M joe.bloggs@sheffield.ac.uk
    #$ -m abe
    #$ -l h_rt=00:10:00
    #$ -l rmem=2G
    #$ -pe smp 4
    module load apps/ansys/19.4/binary
    mapdl -b -np $NSLOTS -smp -i CrankSlot_Flexible.inp

The equivalent batch script for using MPI (``multi-node distributed memory``) parallel environment is ::

    #!/bin/bash
    #$ -cwd
    #$ -N mech_test
    #$ -M joe.bloggs@sheffield.ac.uk
    #$ -m abe
    #$ -l h_rt=00:10:00
    #$ -l rmem=2G
    #$ -pe mpi 4
    module load apps/ansys/19.4/binary
    mapdl -i CrankSlot_Flexible.inp -b -np $NSLOTS -sge -mpi=INTELMPI -rsh -sgepe mpi-rsh 

		
Installation notes
------------------

Ansys 15.0 was installed using the
:download:`install_ansys.sh </sharc/software/install_scripts/apps/ansys/15.0/install_ansys.sh>` script; the module
file is
:download:`/usr/local/modulefiles/apps/ansys/15.0/binary </sharc/software/modulefiles/apps/ansys/15.0/binary>`.

Ansys 16.1 was installed using the
:download:`install_ansys.sh </sharc/software/install_scripts/apps/ansys/16.1/install_ansys.sh>` script; the module
file is
:download:`/usr/local/modulefiles/apps/ansys/16.1 </sharc/software/modulefiles/apps/ansys/16.1>`.

Ansys 17.2 was installed using the
:download:`install_ansys.sh </sharc/software/install_scripts/apps/ansys/17.2/install_ansys.sh>` script; the module
file is
:download:`/usr/local/modulefiles/apps/ansys/17.2 </sharc/software/modulefiles/apps/ansys/17.2>`. 

Ansys 18.0 was installed using the
:download:`install_ansys_180.sh </sharc/software/install_scripts/apps/ansys/18.0/binary/install_ansys_180.sh>` script; the module
file is
:download:`/usr/local/modulefiles/apps/ansys/18.0/binary </sharc/software/modulefiles/apps/ansys/18.0/binary>`. 

Ansys 18.2 was installed using the
:download:`install_ansys_182.sh </sharc/software/install_scripts/apps/ansys/18.2/binary/install_ansys_182.sh>` script; the module
file is
:download:`/usr/local/modulefiles/apps/ansys/18.2/binary </sharc/software/modulefiles/apps/ansys/18.2/binary>`. 

Ansys 19.0 was installed using the
:download:`install_ansys_190.sh </sharc/software/install_scripts/apps/ansys/19.0/binary/install_ansys_190.sh>` script; the module
file is
:download:`/usr/local/modulefiles/apps/ansys/19.0/binary </sharc/software/modulefiles/apps/ansys/19.0/binary>`.

Ansys 19.1 was installed using the
:download:`install_ansys_191.sh </sharc/software/install_scripts/apps/ansys/19.1/binary/install_ansys_191.sh>` script; the module
file is
:download:`/usr/local/modulefiles/apps/ansys/19.1/binary </sharc/software/modulefiles/apps/ansys/19.1/binary>`.

Ansys 19.2 was installed using the
:download:`install_ansys_192.sh </sharc/software/install_scripts/apps/ansys/19.2/binary/install_ansys_192.sh>` script; the module
file is
:download:`/usr/local/modulefiles/apps/ansys/19.2/binary </sharc/software/modulefiles/apps/ansys/19.2/binary>`.

Ansys 19.3 was installed using the
:download:`install_ansys_193.sh </sharc/software/install_scripts/apps/ansys/19.3/binary/install_ansys_193.sh>` script; the module
file is
:download:`/usr/local/modulefiles/apps/ansys/19.3/binary </sharc/software/modulefiles/apps/ansys/19.3/binary>`.

Ansys 19.4 was installed using the
:download:`install_ansys_194.sh </sharc/software/install_scripts/apps/ansys/19.4/binary/install_ansys_194.sh>` script; the module
file is
:download:`/usr/local/modulefiles/apps/ansys/19.4/binary </sharc/software/modulefiles/apps/ansys/19.4/binary>`.

The ``mpi-rsh`` tight-integration parallel environment is required to run Ansys/Fluent using MPI due to 
SSH access to worker nodes being prohibited for most users.

For versions 19.3 & 19.4 mapdl will not run without modifying the file::

    /usr/local/packages/live/eb/ANSYS/19.4/v194/ansys/bin/anssh.ini

The following instruction should be inserted at line 2127 in ``anssh.ini``::

    setenv KMP_AFFINITY compact
