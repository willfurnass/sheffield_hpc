.. _simbody:

SimBody
=======

.. sidebar:: SimBody

   :Version: 3.2.1
   :Support Level: Bronze
   :Dependancies: libs/hdf5/gcc/openmpi/1.8.14
   :URL: http://cgns.github.io/WhatIsCGNS.html
   :Location: /usr/local/packages6/libs/gcc/4.4.7/cgnslib


Usage
-----
To make this library available, run the following module command

.. code-block:: none

        module load libs/gcc/4.4.7/cgns/3.2.1

This will also load the module files for the prerequisite libraries, Open MPI 1.8.3 and HDF5 1.8.14 with parallel support.

Installing
----------
This section is primarily for administrators of the system.

Installed using the following procedure::

    module load compilers/gcc/4.8.2

    cmake ../simbody-Simbody-3.5.3/ -DCMAKE_INSTALL_PREFIX=/usr/local/packages6/libs/gcc/4.4.7/simbody/3.5.3

    make -j 8


Module File
-----------
Module File Location: ``/usr/local/modulefiles/libs/gcc/4.4.7/cgns/3.2.1``

.. code-block:: none

	#%Module1.0#####################################################################
	##
	## cgns 3.2.1 module file
	##

	## Module file logging
	source /usr/local/etc/module_logging.tcl
	##

	proc ModulesHelp { } {
		puts stderr "Makes the cgns 3.2.1 library available"
	}

	module-whatis   "Makes the cgns 3.2.1 library available"
	module load libs/hdf5/gcc/openmpi/1.8.14

	set CGNS_DIR /usr/local/packages6/libs/gcc/4.4.7/cgnslib/3.2.1

	prepend-path LD_LIBRARY_PATH $CGNS_DIR/lib
	prepend-path CPATH $CGNS_DIR/include
