.. Warning:: 
    Iceberg reaches end-of-life on **30th November 2020.**
    If you are running jobs on Iceberg then you need to take urgent action to ensure that your jobs/scripts will run on ShARC or Bessemer. 
 
    If you have never used ShARC or Bessemer then now is the time to test your scripts.
    Not all software on Iceberg is available on ShARC/Bessemer. 

.. _beagle:

beagle
======

.. sidebar:: beagle

   :Latest version: 2.1.2
   :URL: https://github.com/beagle-dev/beagle-lib
   :Location: /usr/local/packages6/libs/gcc/4.4.7/beagle/2.1.2/

General purpose library for evaluating the likelihood of sequence evolution on trees.

Usage
-----
To make this library available, run the following module command: ::

      module load libs/gcc/4.4.7/beagle/2.1.2

This populates the environment variables ``C_INCLUDE_PATH``, ``LD_LIBRARY_PATH`` and ``LD_RUN_PATH`` with the relevant directories.

Installation notes
------------------
This section is primarily for administrators of the system.

Beagle 2.1.2 was compiled with gcc 4.4.7.

* :download:`Install script </iceberg/software/install_scripts/libs/gcc/4.4.7/beagle/2.1.2/install.sh>` 
* :download:`Module file </iceberg/software/modulefiles/libs/gcc/4.4.7/beagle/2.1.2>`
