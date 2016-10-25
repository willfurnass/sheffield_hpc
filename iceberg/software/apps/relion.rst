RELION
======

.. sidebar:: relion

   :Versions:  2.0.1-beta, 1.4
   :URL: http://www2.mrc-lmb.cam.ac.uk/relion/

RELION is a software package that performs an empirical Bayesian approach to (cryo-EM) structure determination by single-particle analysis. Note that RELION is distributed under a GPL license. 

Making RELION available
-----------------------
The following module command makes the latest version of RELION available to your session ::

      module load apps/gcc/4.4.7/relion

Alternatively, you can make a specific version available using one of the following ::

      module load apps/gcc/4.4.7/relion/2.0.1-beta
      module load apps/gcc/4.4.7/relion/1.4

Both versions support parallelisation using multiple cores on single nodes (i.e. a shared memory approach using OpenMP) and cores on different nodes (i.e. a distributed memory approach using MPI) plus the combination of the two approaches.

Version ``2.0.1-beta`` also includes GPU support (using CUDA 7.5).  

Installation notes
------------------
These are primarily for system administrators.

Install instructions: http://www2.mrc-lmb.cam.ac.uk/relion/index.php/Download_%26_install

Version 2.0.1-beta
^^^^^^^^^^^^^^^^^

* Commit ``a6225608e270ac68af463ec3ba7567dbbbc4a4a0`` of `this repository <https://bitbucket.org/tcblab/relion2-beta.git>`_.
* Installed using the GCC 4.4.7 compiler, OpenMPI 1.10.1, CUDA 7.5, FFTW 3.3.4 and FLTK 1.3.3.
    * See `install_relion_2.0.1-beta.sh <https://github.com/mikecroucher/HPC_Installers/blob/master/apps/relion/2.0.1-beta/sheffield/iceberg/install_relion_2.0.1-beta.sh>`_
* `This modulefile <https://github.com/mikecroucher/HPC_Installers/blob/master/apps/relion/2.0.1-beta/sheffield/iceberg/relion_2.0.1-beta_modulefile>`_ was installed as ``/usr/local/modulefiles/apps/gcc/4.4.7/relion/2.0.1-beta``

Version 1.4
^^^^^^^^^^^

* Installed using the GCC 4.4.7 compiler and OpenMPI 1.8.8
* `install_relion.sh <https://github.com/rcgsheffield/sheffield_hpc/blob/master/software/install_scripts/apps/gcc/4.4.7/relion/1.4/install_relion.sh>`_
>>>>>>> Iceberg: start adding info on RELION 2.0.1-beta
* Note - the environment variable RELION_QSUB_TEMPLATE points to an SGE qsub template, which needs customizing to work with our environment

The module file is on the system at `/usr/local/modulefiles/apps/gcc/4.4.7/relion/1.4`

The contents of the module file is ::

  #%Module1.0#####################################################################
  ##
  ## relion module file
  ##

  ## Module file logging
  source /usr/local/etc/module_logging.tcl
  ##

  proc ModulesHelp { } {
      global bedtools-version

      puts stderr "   Setups `relion-$relionversion' environment variables"
  }

  set     relionversion 1.4

  module load apps/gcc/4.4.7/ctffind/3.140609
  module load apps/binapps/resmap/1.1.4
  module load mpi/gcc/openmpi/1.8.8

  prepend-path PATH /usr/local/packages6/apps/gcc/4.4.7/relion/1.4/bin
  prepend-path LD_LIBRARY_PATH /usr/local/packages6/apps/gcc/4.4.7/relion/1.4/lib
  setenv RELION_QSUB_TEMPLATE /usr/local/packages6/apps/gcc/4.4.7/relion/1.4/bin/relion_qsub.csh
  setenv RELION_CTFFIND_EXECUTABLE ctffind3_mp.exe
  setenv RELION_RESMAP_EXECUTABLE /usr/local/packages6/apps/binapps/resmap/1.1.4
