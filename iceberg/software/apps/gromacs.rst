.. Warning:: 
    Iceberg reaches end-of-life on **30th November 2020.**
    If you are running jobs on Iceberg then you need to take urgent action to ensure that your jobs/scripts will run on ShARC or Bessemer. 
 
    If you have never used ShARC or Bessemer then now is the time to test your scripts.
    Not all software on Iceberg is available on ShARC/Bessemer. 

.. _gromacs:

GROMACS
=======

.. sidebar:: GROMACS

   :Latest version: 5.1
   :Dependancies: mpi/intel/openmpi/1.10.0
   :URL: http://www.gromacs.org/


Interactive Usage
-----------------
After connecting to iceberg (see :ref:`ssh`),  start an interactive session with the ``qsh`` or ``qrsh`` command. 
To make GROMACS available in this session, run one of the following command: ::

      source /usr/local/packages6/apps/intel/15/gromacs/5.1/bin/GMXRC

Installation notes
-------------------

Version 5.1
###########

Compilation Choices:

* Use latest intel compilers
* ``-DGMX_MPI=on``
* ``-DCMAKE_PREFIX_PATH=/usr/local/packages6/apps/intel/15/gromcas``
* ``-DGMX_FFT_LIBRARY="fftw3"``
* ``-DGMX_BUILD_OWN_FFTW=ON``

The script used to build gromacs can be found :download:`here
</iceberg/software/install_scripts/apps/gromacs/install_gromacs_5.1.sh>`.
