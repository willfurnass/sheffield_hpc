VASP
====

.. sidebar:: VASP

   :Version: 5.4.4
   :Dependencies: Fortran and C compilers, an implementation of MPI, numerical libraries BLAS, LAPACK, ScaLAPACK, FFTW. Modules for Intel compiler 2022.2.1, Intel MPI 2021.7.1 and Intel MKL 2022.2.1 loaded.
   :URL: https://www.vasp.at/
   :Documentation: https://www.vasp.at/documentation


The Vienna Ab initio Simulation Package (VASP) is a computer program for atomic scale materials modelling, e.g. electronic structure calculations and quantum-mechanical molecular dynamics, from first principles.


Usage
-----

VASP 5.4.4 can be activated using one of the module load commands below: ::

    source /opt/apps/testapps/el7/software/staging/Lmod/7.3/lmod/7.3/init/bash
    module use /opt/apps/testapps/el7/modules/staging/all/
    module load VASP/5.4.4-intel-2022b

The VASP executables are ``vasp_std``, ``vasp_gam`` and ``vasp_ncl``.

.. important::

    Only licensed users of VASP are entitled to use the code; refer to VASP's website for license details: https://www.vasp.at/registration_form/ . Access to VASP on Bessemer is restricted to members of the unix group ``hpc_vasp``.
    To be added to this group, please contact ``research-it@sheffield.ac.uk`` and provide evidence of your eligibility to use VASP.


Pseudopotential files
^^^^^^^^^^^^^^^^^^^^^

VASP Pseudopotentials can be found in the ``/opt/apps/testapps/common/shared/VASP/VASP_POTCAR/`` directory in versioned ``potpaw`` sub-directories after loading the module, e.g. :

.. code-block:: console

    $ ls /opt/apps/testapps/common/shared/VASP/VASP_POTCAR/5.4
    potpaw_LDA  potpaw_PBE

Batch jobs
----------

Users are encouraged to write their own batch submission scripts. The following is an example batch submission script, ``run.sh``, to run ``vasp_std`` and which is submitted to the queue with the command ``sbatch run.sh``. ::

    #!/bin/bash
    #SBATCH --nodes 16
    #SBATCH --ntasks-per-node 1
    #SBATCH --mem=16000
    #SBATCH --time 00:30:00
    #SBATCH --mail-user=jane.doe@sheffield.ac.uk
    #SBATCH --mail-type=ALL
    source /opt/apps/testapps/el7/software/staging/Lmod/7.3/lmod/7.3/init/bash
    module use /opt/apps/testapps/el7/modules/staging/all/
    module load VASP/5.4.4-intel-2022b
    export PSEUDOPOTENTIAL_DIR=/opt/apps/testapps/common/shared/VASP/VASP_POTCAR/5.4

    rm results.dat
    drct=$(pwd)

    for i in std_relaxation constrMD_microcanonical constrMD_canonical
    do
      cd $drct/$i
      ln -sf ../POTCAR .
      ln -sf ../POSCAR .
      ln -sf ../KPOINTS .
      srun --export=ALL vasp_std
    /bin/rm CHG* WAVECAR
    done

The script requests 16 cores (one core per node) using Intel MPI with a runtime of 30 mins and 16 GB of real memory. The batch script above with the input for the "Adsorption of H2O on TiO2" example (https://www.vasp.at/wiki/index.php/Adsorption_of_H2O_on_TiO2) from the online VASP tutorials (https://www.vasp.at/wiki/index.php/Category:Tutorials).


Installation notes
------------------

Not relevant for Pilot User phase.
