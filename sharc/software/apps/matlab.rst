.. _matlab_sharc:

MATLAB
======

.. sidebar:: MATLAB

   :Versions:  2016a 2016b
   :Support Level: FULL
   :Dependancies: None
   :URL: http://uk.mathworks.com/products/matlab
   :Documentation: http://uk.mathworks.com/help/matlab

Scientific computing and visualisation.

Interactive Usage
-----------------
After connecting to ShARC (see :ref:`ssh`),  start an interactive session with the ``qsh`` command.

The latest version of MATLAB (currently 2016b) is made available by running: ::

        module load apps/matlab

Alternatively, you can load a specific version with one of of the following commands: ::

        module load apps/matlab/2016a/binary
        module load apps/matlab/2016b/binary

You can then run MATLAB by entering ``matlab``

Serial (one CPU) Batch usage
----------------------------
Here, we assume that you wish to run the program ``hello.m`` on the system.

First, you need to write a batch submission file. We assume you'll call this ``my_job.sge``: ::

        #!/bin/bash
        #$ -l rmem=4G                        # Request  4 GB of real memory
        #$ -l mem=16G                        # Request 16 GB of virtual memory
        $ -cwd                               # Run job from current directory
        module load apps/matlab/2016b/binary # Make specific version of MATLAB available

        matlab -nodesktop -r 'hello'

Ensuring that ``hello.m`` and ``my_job.sge`` are both in your current working directory, submit your job to the batch system: ::

        qsub my_job.sge

Note that we are running the script ``hello.m`` but we drop the ``.m`` in the call to MATLAB. That is, we do ``-r 'hello'`` rather than ``-r hello.m``.

On the `:ref:iceberg cluster<matlab_iceberg>` MATLAB batch jobs can also be submitted using a script called ``runmatlab``; this script has been deprecated on this cluster.

MATLAB Compiler and running free-standing compiled MATLAB programs
------------------------------------------------------------------

The MATLAB compiler **mcc** can be used to generate standalone executables.
These executables can then be run on other computers that does not have MATLAB installed. 
We strongly recommend you use R2016a or later versions to take advantage of this feature. 

To compile a MATLAB function or script for example called ``myscript.m`` the following steps are required: ::

        # Load the matlab 2016a module
        module load apps/matlab/2016a/binary  

        # Compile your program to generate the executable myscript and 
        # also generate a shell script named run_myscript.sh 
        mcc -m myscript.m

        # Finally run your program
        ./run_myscript.sh $MCRROOT

If ``myscript.m`` is a MATLAB function that require inputs these can be suplied on the command line. 
For example if the first line of ``myscript.m`` reads: ::

        function out = myscript ( a , b , c )

then to run it with 1.0, 2.0, 3.0 as its parameters you will need to type: ::

        ./run_myscript.sh $MCRROOT 1.0 2.0  3.0 

After a successful compilation and running you can transfer your executable and the runscript to another computer.
That computer does not have to have MATLAB installed or licensed on it but it will have to have the MATLAB runtime system installed. 
This can be done by either downloading the MATLAB runtime environment from Mathworks web site or 
by copying the installer file from the cluster itself which resides in: ::

        $MCRROOT/toolbox/compiler/deploy/glnxa64/MCRInstaller.zip

This file can be unzipped in a temporary area and run the setup script that unzipping yields to install the MATLAB runtime environment.
Finally the environment variable ``$MCRROOT`` can be set to the directory containing the runtime environment.  
 
Parallel MATLAB
---------------

**Not yet configured on this cluster.**

Training
--------

* CiCS run an `Introduction to Matlab course <http://rcg.group.shef.ac.uk/courses/matlab/>`_
* In November 2015, CiCS hosted a masterclass in *Parallel Computing in MATLAB*. The materials `are available online <http://rcg.group.shef.ac.uk/courses/mathworks-parallelmatlab/>`_


Installation notes
------------------

These notes are primarily for system administrators.

Installation and configuration is a four-stage process:

* Set up the floating license server (the license server for earlier MATLAB versions can be used), ensuring that it can serve licenses for any new versions of MATLAB that you want to install
* Run a graphical installer to download MATLAB *archive* files used by the main (automated) installation process
* Run the same installer in 'silent' command-line mode to perform the installation using those archive files and a text config file.
* Install a relevant modulefile

In more detail:

#. If necessary, update the floating license keys on ``licserv4.shef.ac.uk`` to ensure that the licenses are served for the versions to install.
#. Log on to Mathworks site to download the MATLAB installer package for 64-bit Linux ( for R2016a this was called ``matlab_R2016b_glnxa64.zip`` )

#. ``unzip`` the installer package in a directory with ~10GB of space (needed as many MATLAB *archive* files will subsequently be downloaded here).  Using a directory on an NFS mount (e.g. ``/data/${USER}/MathWorks/R2016b``) allows the same downloaded archives to be used to install MATLAB on multiple clusters.
#. ``./install`` to start the graphical installer (needed to download the MATLAB archive files).
#. Select install choice of *Log in to Mathworks Account* and log in with a *License Administrator* account (not a *Licensed End User* (personal) account).
#. Select *Download only*.
#. Select the offered default *Download path* and select the directory you ran ``./install`` from.  Wait a while for all requested archive files to be downloaded.
#. Next, ensure ``installer_input.txt`` looks like the following ::
    
    fileInstallationKey=XXXXX-XXXXX-XXXXX-XXXXX-XXXXX-XXXXX-XXXXX-XXXXX-XXXXX-XXXXX-XXXXX-XXXXX-XXXXX-XXXXX-XXXXX-XXXXX-XXXXX-XXXXX-XXXXX-XXXXX-XXXXX
    agreeToLicense=yes
    outputFile=matlab_2016b_install.log
    mode=silent
    licensePath=/usr/local/packages/matlab/network.lic
    lmgrFiles=false
    lmgrService=false

#. Create the installation directory e.g.: ::

    mkdir -m 2755 -p /usr/local/packages/apps/matlab/R2016b/binary
    chown ${USER}:app-admins /usr/local/packages/apps/matlab/R2016b/binary

#. Run the installer using our customized ``installer_input.txt`` like so: ``./install -mode silent -inputFile ${PWD}/installer_input.txt`` ; installation should finish with exit status ``0`` if all has worked.
#. Install a *modulefile* to prepend to to the ``PATH`` environment variable and set the ``MCRROOT`` environment variable (used by the ``mcc`` compiler):
    
    - :download:`This modulefile </sharc/software/modulefiles/apps/matlab/2016b/binary>` was installed as ``/usr/local/modulefiles/apps/matlab/2016b/binary``
    - :download:`This modulefile </sharc/software/modulefiles/apps/matlab/2016a/binary>` was installed as ``/usr/local/modulefiles/apps/matlab/2016a/binary``

#. Ensure the contents of the install directory and the modulefile are writable by those in ``app-admins`` group e.g.: ::

    chmod -R g+w ${USER}:app-admins /usr/local/packages6/matlab/R2016b /usr//local/modulefiles/apps/matlab/2016b

**TODO**: Documentation for MATLAB parallel configuration.
