HDF5
====

.. sidebar:: HDF5
   
   :Maintainer: Stuart Mumford
   :Email: stuart.mumford@sheffield.ac.uk
   :Version: 1.8.13
   :Support Level: bronze
   :Dependancies: openmpi (1.8.13)
   :URL: http://www.hdfgroup.org/HDF5/
   :Documentation: http://www.hdfgroup.org/HDF5/doc/ 


HDF5 is a data model, library, and file format for storing and managing data.
It supports an unlimited variety of datatypes, and is designed for flexible and efficient I/O and for high volume and complex data.
HDF5 is portable and is extensible, allowing applications to evolve in their use of HDF5.
The HDF5 Technology suite includes tools and applications for managing, manipulating, viewing, and analyzing data in the HDF5 format. 

Usage
-----

A module file is provided along with this library, it can be enabled using the standard module loading syntax::

     module load libs/hdf5/gcc/openmpi/1.8.13
     module load libs/hdf5/gcc/1.8.13


Installing
----------

This package is built from the source code distribution from the HDF Group website.

Two primary versions of this library are provided, a MPI parallel enabled version and a serial version.
The serial version has the following configuration flags enabled::

    --enable-fortran --enable-fortran2003 --enable-cxx --enable-shared

The parallel version has the following flags::

    --enable-fortran --enable-fortran2003 --enable-shared --enable-parallel

The parallel library does not support C++, hence it being disabled for the parallel build.
