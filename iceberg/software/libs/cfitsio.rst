.. Warning:: 
    Iceberg reaches end-of-life on **30th November 2020.**
    If you are running jobs on Iceberg then you need to take urgent action to ensure that your jobs/scripts will run on ShARC or Bessemer. 
 
    If you have never used ShARC or Bessemer then now is the time to test your scripts.
    Not all software on Iceberg is available on ShARC/Bessemer. 

.. _cfitsio:

cfitsio
=======

.. sidebar:: cfitsio

   :Latest version: 3.380
   :URL: http://heasarc.gsfc.nasa.gov/fitsio/fitsio.html

CFITSIO is a library of C and Fortran subroutines for reading and writing data
files in FITS (Flexible Image Transport System) data format. CFITSIO provides
simple high-level routines for reading and writing FITS files that insulate
the programmer from the internal complexities of the FITS format. 

Usage
-----
To make this library available, run the following module command ::

        module load libs/gcc/5.2/cfitsio

The modulefile creates a variable ``$CFITSIO_INCLUDE_PATH`` which is the path
to the include directory.

Installation notes
------------------
This section is primarily for administrators of the system. CFITSIO 3.380 was compiled with gcc 5.2.
The compilation used this :download:`script </iceberg/software/install_scripts/libs/gcc/5.2/install_cfitsio.sh>` 
and it is loaded with this :download:`modulefile </iceberg/software/modulefiles/libs/gcc/5.2/cfitsio/3.380>`.
