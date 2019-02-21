GROMACS
=======

.. sidebar:: GROMACS

   :Version: 2016.4, 2018.1
   :Dependencies: Modules loaded GCC 4.9.4, CUDA 8.0.44 (GPU support), and OPENMPI 2.0.1 (MPI support).
   :URL: http://www.gromacs.org/
   :Documentation: http://manual.gromacs.org/documentation/
   :Tutorials: http://www.gromacs.org/Documentation/Tutorials


GROMACS is a versatile package to perform molecular dynamics, i.e. simulate the Newtonian equations of motion for systems with hundreds to millions of particles.
It is primarily designed for biochemical molecules like proteins, lipids and nucleic acids that have a lot of complicated bonded interactions, but since GROMACS is extremely fast at calculating the nonbonded interactions (that usually dominate simulations) many groups are also using it for research on non-biological systems, e.g. polymers.


Usage
-----

GROMACS 2016.4 can be activated using the module files::

    module load apps/gromacs/2016.4/gcc-4.9.4
    module load apps/gromacs/2018.1/gcc-4.9.4
    module load apps/gromacs/2016.4/gcc-4.9.4-cuda-8.0
    module load apps/gromacs/2018.1/gcc-4.9.4-cuda-8.0
    module load apps/gromacs/2018.1/gcc-4.9.4-openmpi-2.0.1

Where the module files ``apps/gromacs/201*/gcc-4.9.4-cuda-8.0`` is for the installation compiled to run on GPUs.
The module file ``apps/gromacs/201*/gcc-4.9.4-openmpi-2.0.1`` is for the installation compiled to run using OPENMPI parallelism.
The GROMACS executable is ``gmx``. Typing ``gmx help commands`` will display a list of commands for ``gmx`` and their function.


Batch jobs
----------

Users are encouraged to write their own batch submission scripts. The following is an example batch submission script, ``my_job.sh``, to run ``gmx`` and which is submitted to the queue by typing ``qsub my_job.sh``. ::

    #!/bin/bash
    #$ -cwd
    #$ -l h_rt=06:00:00
    #$ -l rmem=2G
    #$ -pe smp 1
    #$ -l gpu=1

    module load apps/gromacs/2016.4/gcc-4.9.4-cuda-8.0

    gmx grompp -f grompp.mdp -c conf.gro -p topol.top -o topol.tpr
    gmx mdrun -s topol.tpr

The script requests one CPU core using the OpenMP parallel environment ``smp``, with 2 GB of real memory per CPU core, and one GPU per one CPU core. The requested runtime is 6 hours.
The GROMACS input line beginning with ``gmx grompp`` is to make a run input file; ``gmx mdrun`` is to perform a simulation, do a normal mode analysis or an energy minimization.


Installation notes
------------------

Four GROMACS 2016.4 & 2018.1 installations are available on ShARC; two with and two without GPU support. Both installations use single-node OpenMP parallelism, are single-precision builds and use the GROMACS installation of the FFTW3 library. One installation supports OPENMPI parallelism (note: request -pe mpi no_of_cores in your batch script rather than -pe smp no_of_cores).

GROMACS 2016.4 without GPU support was installed using the
:download:`install_gromacs.sh </sharc/software/install_scripts/apps/gromacs/2016.4/gcc-4.9.4/install_gromacs.sh>` script;
the module file is
:download:`/usr/local/modulefiles/apps/gromacs/2016.4/gcc-4.9.4 </sharc/software/modulefiles/apps/gromacs/2016.4/gcc-4.9.4>`.

GROMACS 2018.1 without GPU support was installed using the
:download:`install_gromacs.sh </sharc/software/install_scripts/apps/gromacs/2018.1/gcc-4.9.4/install_gromacs.sh>` script;
the module file is
:download:`/usr/local/modulefiles/apps/gromacs/2018.1/gcc-4.9.4 </sharc/software/modulefiles/apps/gromacs/2018.1/gcc-4.9.4>`.

GROMACS 2016.4 with GPU support was installed using the
:download:`install_gromacs_gpu.sh </sharc/software/install_scripts/apps/gromacs/2016.4/gcc-4.9.4-cuda-8.0/install_gromacs_gpu.sh>` script;
the module file is
:download:`/usr/local/modulefiles/apps/gromacs/2016.4/gcc-4.9.4-cuda-8.0 </sharc/software/modulefiles/apps/gromacs/2016.4/gcc-4.9.4-cuda-8.0>`.

GROMACS 2018.1 with GPU support was installed using the
:download:`install_gromacs_gpu.sh </sharc/software/install_scripts/apps/gromacs/2018.1/gcc-4.9.4-cuda-8.0/install_gromacs_gpu.sh>` script;
the module file is
:download:`/usr/local/modulefiles/apps/gromacs/2018.1/gcc-4.9.4-cuda-8.0 </sharc/software/modulefiles/apps/gromacs/2018.1/gcc-4.9.4-cuda-8.0>`.

GROMACS 2018.1 with OPENMPI support was installed using the
:download:`install_gromacs_gpu.sh </sharc/software/install_scripts/apps/gromacs/2018.1/gcc-4.9.4-openmpi-2.0.1/install_gromacs_mpi.sh>` script;
the module file is
:download:`/usr/local/modulefiles/apps/gromacs/2018.1/gcc-4.9.4-openmpi-2.0.1 </sharc/software/modulefiles/apps/gromacs/2018.1/gcc-4.9.4-openmpi-2.0.1>`.
The GROMACS 2016.4 & 2018.1 installations were tested by using ``make check`` to run regression tests as part of the installation process.
