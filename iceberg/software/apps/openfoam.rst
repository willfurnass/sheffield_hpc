.. include:: warning.rst 


OpenFOAM
========

.. sidebar:: OpenFOAM
   
   :Versions: 3.0.0, 2.4.0
   :Dependencies: gcc/5.2 mpi/gcc/openmpi/1.10.0 libs/gcc/5.2/scotch/6.0.4
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

These ``module load`` commands also load:

* :ref:`OpenMPI <openmpi_gcc_iceberg>` 1.10.0 
* :ref:`GCC <gcc_iceberg>` 5.2
* the :ref:`Scotch <scotch>` 6.0.4 library

Installation notes
------------------

Version 3.0.0
^^^^^^^^^^^^^

Compiled using this :download:`install_openfoam.sh</iceberg/software/install_scripts/apps/gcc/5.2/openfoam/install_openfoam.sh>` script.

The :download:`3.0.0 module file</iceberg/software/modulefiles/apps/gcc/5.2/openfoam/3.0.0>`.

Version 2.4.0
^^^^^^^^^^^^^

The compilation process was not documented.

The :download:`2.4.0 module file</iceberg/software/modulefiles/apps/gcc/5.2/openfoam/2.4.0>`.
