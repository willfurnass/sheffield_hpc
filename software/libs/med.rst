.. _MED:

MED
===

.. sidebar:: MED

   :Version: 3.0.8
   :Support Level: Bronze
   :Dependancies: None
   :URL: http://www.salome-platform.org/downloads/current-version
   :Location: /usr/local/packages6/libs/gcc/4.4.7/med/3.0.8

The purpose of the MED module is to provide a standard for storing and recovering computer data associated to numerical meshes and fields, and to facilitate the exchange between codes and solvers. 

Usage
-----
To make this library available, run the following module command

.. code-block:: none

        module load libs/gcc/4.4.7/med/3.0.8

Installing
----------
* This is a pre-requisite for Code Saturne version 4.0.
* It was built with gcc 4.4.7

.. code-block:: none

        module load mpi/gcc/openmpi/1.8.3
	tar -xvzf med-3.0.8.tar.gz
	cd med-3.0.8
	mkdir -p /usr/local/packages6/libs/gcc/4.4.7/med/3.0.8
	./configure --prefix=/usr/local/packages6/libs/gcc/4.4.7/med/3.0.8 --disable-fortran --with-hdf5=/usr/local/packages6/hdf5/gcc-4.4.7/openmpi-1.8.3/hdf5-1.8.14/ --disable-python
        make
	make install

Fortran was disabled because otherwise the build failed with compilation errors. It's not needed for Code Saturne 4.0.

Python was disabled because it didn't have MPI support.

testing
-------
.. code-block:: none

        make check

All tests passed
