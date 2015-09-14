.. _hdf5pgi:

HDF5 (PGI build)
================

.. sidebar:: HDF5 (PGI build)

   :Latest version: 1.8.15-patch1
   :Dependancies: PGI Compiler 15.7, PGI Openmpi 1.8.8
   :URL: http://www.hdfgroup.org/HDF5/
   :Documentation: http://www.hdfgroup.org/HDF5/doc/


HDF5 is a data model, library, and file format for storing and managing data.
It supports an unlimited variety of datatypes, and is designed for flexible and efficient I/O and for high volume and complex data.
HDF5 is portable and is extensible, allowing applications to evolve in their use of HDF5.
The HDF5 Technology suite includes tools and applications for managing, manipulating, viewing, and analyzing data in the HDF5 format.

Usage
-----
This version of HDF5 was compiled using version 15.7 of the PGI Compiler and OpenMPI 1.8.8. To make it available, run the following module command after starting a ``qsh`` or ``qrsh`` session ::

    module load libs/hdf5/pgi/1.8.15-patch1

This module also loads the relevant modules for OpenMPI and PGI Compiler. To see which modules have been loaded, use the command ``module list``

Installation notes
------------------
This section is primarily for administrators of the system.

**Version 1.8.15-patch1**

Compiled using PGI 15.7 and OpenMPI 1.8.8

* `install_pgi_hdf5_1.8.15-patch1.sh   <https://github.com/rcgsheffield/blob/master/software/install_scripts/libs/pgi/hdf5/install_pgi_hdf5_1.8.15-patch1.sh>`_ Install script
* `Modulefile <https://github.com/cgsheffield/iceberg_software/blob/master/software/modulefiles/libs/pgi/hdf5/1.8.15-patch1>`_ located on the system at ``/usr/local/modulefiles/libs/hdf5/pgi/1.8.15-patch1``
