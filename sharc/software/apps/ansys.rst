ANSYS
=====

.. sidebar:: ANSYS

   :Versions: 15.0, 16.1, 17.2, 18.0, 18.2, 19.0, 19.1, 19.2, 19.3, 19.4, 20.1, 20.2 & 21.1
   :Dependencies: No prerequsite modules loaded. However, if using the User Defined Functions (UDF) will also need the following: For ANSYS Mechanical, Workbench, CFX and AutoDYN: Intel 14.0 or above; Compiler For Fluent: GCC 4.6.1 or above
   :URL: http://www.ansys.com


The ANSYS suite of programs can be used to numerically simulate a large variety of structural and fluid dynamics problems found in many engineering, physics, medical, aeronotics and automative industry applications.


Usage
-----

ANSYS can be activated using the following module files::

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
    module load apps/ansys/20.1/binary
    module load apps/ansys/20.2/binary
    module load apps/ansys/21.1/binary


The ANSYS Workbench GUI executable is ``ansyswb``. ``ansyswb`` can be launched during an interactive session with X Window support (e.g. an interactive ``qrshx`` session).
The ANSYS Mechanical executable is ``ansys-mechanical`` and ``fluent`` is the executable for Fluent.

------------

NOTE: for ANSYS versions >= 18.0 using the command ``fluent`` results in an unresponsive fluent launcher. To launch fluent and bypass the launcher use ``fluent dim`` where dim = 2d, 2ddp, 3d or 3ddp or unset the following environment variables before running the command::

    unset SGE_TASK_ID
    unset RESTARTED

ANSYS example models
--------------------

ANSYS contains a large number of example models which can be used to become familiar with the software.
The models can be found in::

    /usr/local/packages/apps/ansys/21.1/binary/v211/ansys/data
    /usr/local/packages/apps/ansys/20.2/binary/v202/ansys/data
    /usr/local/packages/apps/ansys/20.1/binary/v201/ansys/data
    /usr/local/packages/apps/ansys/19.4/binary/v194/ansys/data
    /usr/local/packages/apps/ansys/19.3/binary/v193/ansys/data
    /usr/local/packages/apps/ansys/19.2/binary/v192/ansys/data
    /usr/local/packages/apps/ansys/19.1/binary/v191/ansys/data
    /usr/local/packages/apps/ansys/19.0/binary/v190/ansys/data
    /usr/local/packages/apps/ansys/18.2/binary/v182/ansys/data
    /usr/local/packages/apps/ansys/18.0/binary/v180/ansys/data
    /usr/local/packages/apps/ansys/17.2/v172/ansys/data
    /usr/local/packages/apps/ansys/16.1/v161/ansys/data
    /usr/local/packages/apps/ansys/15.0.7/ansys_inc/v150/ansys/data

Batch jobs
----------
ANSYS Fluent
#############
``Fluent CFD``: the following is an example batch submission script, ``cfd_job.sh``, to run the executable ``fluent`` with input journal file ``test.jou``, and carry out a 2D double precision CFD simulation. The script requests 8 cores using the MPI parallel environment ``mpi-rsh`` with a runtime of 30 mins and 2 GB of real memory per core. The Fluent input journal file is ``test.jou``. **Note:** Use of the ``mpi-rsh`` parallel environment to run MPI parallel jobs for Ansys is required. ::

    #!/bin/bash
    #$ -V
    #$ -cwd
    #$ -M joe.bloggs@sheffield.ac.uk
    #$ -m abe
    #$ -l h_rt=00:30:00
    #$ -l rmem=2G
    #$ -pe mpi-rsh 8
    #$ -N JobName

    module load apps/ansys/20.2/binary

    fluent 2ddp -i test.jou -g -t$NSLOTS -mpi=intel -rsh -sgepe mpi-rsh -sge -driver null

.. note::

    **$NSLOTS** is a Sun Grid Engine variable which will return the requested number of cores, **-rsh** tells Fluent to use RSH instead of SSH, **-sge** forces Fluent to recognise job submission via SGE, **-sgepe** selects the *mpi-rsh* SGE parallel environment and **-driver null** instructs Fluent that it will be running with no GUI to avoid errors caused by plot / figure export.

The job is submitted to the queue by typing::

    qsub cfd_job.sh

|

------------

ANSYS Mechnical / Map-DL
#########################
``Mapdl mechanical``: the following is an example batch submission script, ``mech_job.sh``, to run the mechanical executable ``mapdl`` with input file ``CrankSlot_Flexible.inp``, and carry out a mechanical simulation. The script requests 4 cores using the OpenMP (``single node shared memory``) parallel environment with a runtime of 10 mins and 2 GB of real memory per core. ::

    #!/bin/bash
    #$ -V
    #$ -cwd
    #$ -N JobName
    #$ -M joe.bloggs@sheffield.ac.uk
    #$ -m abe
    #$ -l h_rt=00:10:00
    #$ -l rmem=2G
    #$ -pe smp 4
    module load apps/ansys/20.2/binary
    mapdl -b -np $NSLOTS -smp -i CrankSlot_Flexible.inp

**Note:** Use of the ``mpi-rsh`` parallel environment to run MPI parallel jobs for Ansys is required.
The equivalent batch script for using MPI (``multi-node distributed memory``) parallel environment is ::

    #!/bin/bash
    #$ -V
    #$ -cwd
    #$ -N JobName
    #$ -M joe.bloggs@sheffield.ac.uk
    #$ -m abe
    #$ -l h_rt=00:10:00
    #$ -l rmem=2G
    #$ -pe mpi-rsh 4
    module load apps/ansys/20.2/binary
    mapdl -i CrankSlot_Flexible.inp -b -np $NSLOTS -sge -mpi=INTELMPI -rsh -sgepe mpi-rsh

.. note::

        **$NSLOTS** is a Sun Grid Engine variable which will return the requested number of cores, **-rsh** tells Mechanical to use RSH instead of SSH, **-sge** forces Mechanical to recognise job submission via SGE and **-sgepe** selects the *mpi-rsh* SGE parallel environment.


Installation notes
------------------

ANSYS 15.0 was installed using the
:download:`install_ansys.sh </sharc/software/install_scripts/apps/ansys/15.0/install_ansys.sh>` script; the module
file is
:download:`/usr/local/modulefiles/apps/ansys/15.0/binary </sharc/software/modulefiles/apps/ansys/15.0/binary>`.

ANSYS 16.1 was installed using the
:download:`install_ansys.sh </sharc/software/install_scripts/apps/ansys/16.1/install_ansys.sh>` script; the module
file is
:download:`/usr/local/modulefiles/apps/ansys/16.1 </sharc/software/modulefiles/apps/ansys/16.1>`.

ANSYS 17.2 was installed using the
:download:`install_ansys.sh </sharc/software/install_scripts/apps/ansys/17.2/install_ansys.sh>` script; the module
file is
:download:`/usr/local/modulefiles/apps/ansys/17.2 </sharc/software/modulefiles/apps/ansys/17.2>`.

ANSYS 18.0 was installed using the
:download:`install_ansys_180.sh </sharc/software/install_scripts/apps/ansys/18.0/binary/install_ansys_180.sh>` script; the module
file is
:download:`/usr/local/modulefiles/apps/ansys/18.0/binary </sharc/software/modulefiles/apps/ansys/18.0/binary>`.

ANSYS 18.2 was installed using the
:download:`install_ansys_182.sh </sharc/software/install_scripts/apps/ansys/18.2/binary/install_ansys_182.sh>` script; the module
file is
:download:`/usr/local/modulefiles/apps/ansys/18.2/binary </sharc/software/modulefiles/apps/ansys/18.2/binary>`.

ANSYS 19.0 was installed using the
:download:`install_ansys_190.sh </sharc/software/install_scripts/apps/ansys/19.0/binary/install_ansys_190.sh>` script; the module
file is
:download:`/usr/local/modulefiles/apps/ansys/19.0/binary </sharc/software/modulefiles/apps/ansys/19.0/binary>`.

ANSYS 19.1 was installed using the
:download:`install_ansys_191.sh </sharc/software/install_scripts/apps/ansys/19.1/binary/install_ansys_191.sh>` script; the module
file is
:download:`/usr/local/modulefiles/apps/ansys/19.1/binary </sharc/software/modulefiles/apps/ansys/19.1/binary>`.

ANSYS 19.2 was installed using the
:download:`install_ansys_192.sh </sharc/software/install_scripts/apps/ansys/19.2/binary/install_ansys_192.sh>` script; the module
file is
:download:`/usr/local/modulefiles/apps/ansys/19.2/binary </sharc/software/modulefiles/apps/ansys/19.2/binary>`.

ANSYS 19.3 was installed using the
:download:`install_ansys_193.sh </sharc/software/install_scripts/apps/ansys/19.3/binary/install_ansys_193.sh>` script; the module
file is
:download:`/usr/local/modulefiles/apps/ansys/19.3/binary </sharc/software/modulefiles/apps/ansys/19.3/binary>`.

ANSYS 19.4 was installed using the
:download:`install_ansys_194.sh </sharc/software/install_scripts/apps/ansys/19.4/binary/install_ansys_194.sh>` script; the module
file is
:download:`/usr/local/modulefiles/apps/ansys/19.4/binary </sharc/software/modulefiles/apps/ansys/19.4/binary>`.

----------

ANSYS 20.1, 20.2 and 21.1 were installed using the GUI installer and then permissions were corrected as follows::

    chmod 775 -R /usr/local/packages/apps/ansys/20.1/binary
    chmod 775 -R /usr/local/packages/apps/ansys/20.2/binary
    chmod 775 -R /usr/local/packages/apps/ansys/21.1/binary

Please follow the same install directory structure.

The ``mpi-rsh`` tight-integration parallel environment is required to run ANSYS/Fluent using MPI due to
SSH access to worker nodes being prohibited for most users.

For versions 19.3 & 19.4 and onward mapdl will not run without modifying the file::

    /usr/local/packages/apps/ansys/19.4/binary/v194/ansys/bin/anssh.ini

The following instruction should be inserted at line 2127 in ``anssh.ini``::

    setenv KMP_AFFINITY compact
