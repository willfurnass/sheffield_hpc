
OpenFOAM
========

.. sidebar:: OpenFOAM
   
   :Version: 3.0.0
   :Dependancies: gcc/5.2 mpi/gcc/openmpi/1.10.0 libs/gcc/5.2/scotch/6.0.4
   :URL: http://www.openfoam.org
   :Documentation: http://www.openfoam.org/docs/


OpenFOAM is free, open source software for computational fluid dynamics (CFD).

Usage
-----

You can load a specific version of OpenFOAM using one of: ::

       module load apps/gcc/5.2/openfoam/3.0.0
       module load apps/gcc/5.2/openfoam/2.4.0

This has the same effect as sourcing the OpenFOAM ``bashrc`` file, so that should
not be needed.

These ``module load`` commands also load :ref:`OpenMPI <openmpi_gcc_iceberg>` 1.10.0 and :ref:`GCC <gcc_iceberg>` 5.2.

Next, make sure you also load the :ref:`Scotch <scotch>` library: ::

       module load libs/gcc/5.2/scotch/6.0.4

Installation notes
------------------

OpenFOAM was compiled using the
:download:`install_openfoam.sh</iceberg/software/install_scripts/apps/gcc/5.2/openfoam/install_openfoam.sh>` script.

The :download:`3.0.0 module file</iceberg/software/install_scripts/apps/gcc/5.2/openfoam/3.0.0>`.
