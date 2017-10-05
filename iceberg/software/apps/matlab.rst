.. _matlab_iceberg:

MATLAB
======

.. sidebar:: MATLAB

   :Versions:  2013a , 2013b , 2014a, 2015a, 2016a
   :Support Level: FULL
   :Dependancies: None
   :URL: http://uk.mathworks.com/products/matlab
   :Local URL:  http://www.shef.ac.uk/wrgrid/software/matlab
   :Documentation: http://uk.mathworks.com/help/matlab

Scientific computing and visualisation.

Interactive Usage
-----------------
After connecting to iceberg (see :ref:`ssh`),  start an interactive session with the ``qsh`` command.

The latest version of MATLAB (currently 2016a) is made available with the command: ::

        module load apps/matlab

Alternatively, you can load a specific version with one of of the following commands: ::

        module load apps/matlab/2013a
        module load apps/matlab/2013b
        module load apps/matlab/2014a
        module load apps/matlab/2015a
        module load apps/matlab/2016a

You can then run MATLAB by entering ``matlab``

Serial (one CPU) Batch usage
----------------------------
Here, we assume that you wish to run the program ``helloworld.m`` on the system: ::
	
	function helloworld
		disp('Hello World!')
	end	

First, you need to write a batch submission file. We assume you'll call this ``my_job.sge``: ::

        #!/bin/bash
        #$ -l rmem=4G                  # Request  4 GB of real memory
        #$ -cwd                        # Run job from current directory
        module load apps/matlab/2016a  # Make specific version of MATLAB available

        matlab -nodesktop -nosplash -r helloworld

Ensuring that ``helloworld.m`` and ``my_job.sge`` are both in your current working directory, submit your job to the batch system: ::

        qsub my_job.sge

Note that we are running the script ``helloworld.m`` but we drop the ``.m`` in the call to MATLAB. That is, we do ``-r helloworld`` rather than ``-r helloworld.m``. The output will be written to the job ``.o`` file when the job finishes.

MATLAB Compiler and running free-standing compiled MATLAB programs
------------------------------------------------------------------

The MATLAB compiler **mcc** can be used to generate standalone executables.
These executables can then be run on other computers that does not have MATLAB installed. 
We strongly recommend you use R2016a or later versions to take advantage of this feature. 

To compile a MATLAB function or script for example called ``myscript.m`` the following steps are required: ::

        # Load the matlab 2016a module
        module load apps/matlab/2016a  

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
 

Parallel MATLAB on iceberg
--------------------------

Currently we recommend the 2015a version of MATLAB for parallel work, and task arrays requiring more than a few hours runtime.

The default cluster configuration named **local** provides parallel working environment by 
using the CPUs of the worker node that is running the current MATLAB session.
Each iceberg worker node can run multiple users' jobs simultaneously. 
Therefore depending on who else is using that node at the time, 
parallel MATLAB jobs can create contentions between jobs and slow them considerably. 
It is therefore advisable to start parallel MATLAB jobs that will use the **local** profile from a parallel SGE job.
For example, to use the local profile with 5 workers, do the following;

Start a parallel OpenMP job with 6 workers: ::

        qsh -pe openmp 6

Run MATLAB in that session and select 5 workers: ::

        matlab
        parpool ('local' , 5 )

The above example will use 5 MATLAB workers on a single iceberg node to run a parallel task.  Note that being granted a multi-core interactive session by the scheduler is dependent on Iceberg's loading at the time of request. It is not guaranteed. However the user can reduce the number of cores requested, which will improve their chances of being granted a multi-core session.

To take advantage of the multiple iceberg nodes, you will need to make use of a parallel cluster profile named ``sge``.
This can be done by issuing a locally provided MATLAB command named ``iceberg`` that imports the
parallel cluster profile named ``sge`` that can take advantage of the SGE scheduler to run
larger parallel jobs.

When using the ``sge`` profile, 
MATLAB will be able to submit multiple MATLAB jobs the the SGE scheduler from within MATLAB itself.  
However, each job will have the default resource requirements unless the following trick is deployed.
For example, during your MATLAB session type: ::

    global sge_params
    sge_params='-l rmem=8G -l h_rt=36:00:00'

to make sure that all the MATLAB batch jobs will use up to 8 GBytes of memory and will not be killed
unless they exceed 36 hours of runtime.


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
#. Log on to Mathworks site to download the MATLAB installer package for 64-bit Linux ( for R2016a this was called ``matlab_R2016a_glnxa64.zip`` )

#. ``unzip`` the installer package in a directory with ~10GB of space (needed as many MATLAB *archive* files will subsequently be downloaded here).  Using a directory on an NFS mount (e.g. ``/data/${USER}/MathWorks/R2016a``) allows the same downloaded archives to be used to install MATLAB on multiple clusters.
#. ``./install`` to start the graphical installer (needed to download the MATLAB archive files).
#. Select install choice of *Log in to Mathworks Account* and log in with a *License Administrator* account (not a *Licensed End User* (personal) account).
#. Select *Download only*.
#. Select the offered default *Download path* and select the directory you ran ``./install`` from.  Wait a while for all requested archive files to be downloaded.
#. Next, ensure ``installer_input.txt`` looks like the following ::
    
    fileInstallationKey=XXXXX-XXXXX-XXXXX-XXXXX-XXXXX-XXXXX-XXXXX-XXXXX-XXXXX-XXXXX-XXXXX-XXXXX-XXXXX-XXXXX-XXXXX-XXXXX-XXXXX-XXXXX-XXXXX-XXXXX-XXXXX
    agreeToLicense=yes
    outputFile=matlab_2016a_install.log
    mode=silent
    licensePath=/usr/local/packages6/matlab/network.lic
    lmgrFiles=false
    lmgrService=false

#. Create the installation directory e.g.: ::

    mkdir -m 2755 -p /usr/local/packages6/matlab/R2016a
    chown ${USER}:app-admins /usr/local/packages6/matlab/R2016a

#. Run the installer using our customized ``installer_input.txt`` like so: ``./install -mode silent -inputFile ${PWD}/installer_input.txt`` ; installation should finish with exit status ``0`` if all has worked.
#. Install a *modulefile* with a name and path like ``/usr//local/modulefiles/apps/matlab/2016a`` and contents like ::

    #%Module1.0#####################################################################

    ## Module file logging
    source /usr/local/etc/module_logging.tcl

    proc ModulesHelp { } {
        global version
        puts stderr "	Makes MATLAB 2016a available for use"
    }
    module-whatis   "Makes MATLAB 2016a available"

    # Do not use other versions at the same time.
    conflict apps/matlab

    set     version        2016a
    set     matlabroot     /usr/local/packages6/matlab/R2016a
    set     mcrroot        /usr/local/packages6/matlab/runtime/R2016a/v901
    prepend-path PATH $matlabroot/bin 
    setenv MCRROOT $mcrroot

#. Ensure the contents of the install directory and the modulefile are writable by those in ``app-admins`` group e.g.: ::

    chmod -R g+w ${USER}:app-admins /usr/local/packages6/matlab/R2016a /usr//local/modulefiles/apps/matlab/2016a

**TODO**: Documentation for MATLAB parallel configuration.
