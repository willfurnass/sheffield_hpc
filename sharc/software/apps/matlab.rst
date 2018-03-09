.. _matlab_sharc:

MATLAB
======

.. sidebar:: MATLAB

   :Versions:  2016a, 2016b, 2017a, 2017b
   :Support Level: FULL
   :Dependancies: None
   :URL: http://uk.mathworks.com/products/matlab
   :Documentation: http://uk.mathworks.com/help/matlab

Scientific computing and visualisation.

Interactive Usage
-----------------
After connecting to ShARC (see :ref:`ssh`),  start an interactive session with the ``qrshx`` command.

The latest version of MATLAB (currently 2017b) is made available by running: ::

	module load apps/matlab

Alternatively, you can load a specific version with one of the following commands: ::

	module load apps/matlab/2016a/binary
	module load apps/matlab/2016b/binary
	module load apps/matlab/2017a/binary
	module load apps/matlab/2017b/binary

You can then run MATLAB by entering ``matlab``

Serial (one CPU) Batch usage
----------------------------
Here, we assume that you wish to run the program ``helloworld.m`` on the system: ::
	
	function helloworld
		disp('Hello World!')
	end	

First, you need to write a batch submission file. We assume you'll call this ``my_job.sge``: ::

	#!/bin/bash
	#$ -l rmem=4G                  		# Request  4 GB of real memory
	#$ -cwd                        		# Run job from current directory
	module load apps/matlab/2017b/binary  	# Make specific version of MATLAB available
  
	matlab -nodesktop -nosplash -r helloworld

Ensuring that ``helloworld.m`` and ``my_job.sge`` are both in your current working directory, submit your job to the batch system: ::

	qsub my_job.sge

Note that we are running the script ``helloworld.m`` but we drop the ``.m`` in the call to MATLAB. That is, we do ``-r helloworld`` rather than ``-r helloworld.m``. The output will be written to the job ``.o`` file when the job finishes.

MATLAB Compiler and running free-standing compiled MATLAB programs
------------------------------------------------------------------

The MATLAB compiler **mcc** can be used to generate standalone executables.
These executables can then be run on other computers that does not have MATLAB installed. 
We strongly recommend you use R2016b or later versions to take advantage of this feature. 

To compile a MATLAB function or script for example called ``myscript.m`` the following steps are required: ::

	# Load the matlab 2017b module
	module load apps/matlab/2017b/binary  

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
 
Parallel MATLAB - Single node
-----------------------------

Parallel Matlab can be run exclusively on a single node (using a maximum of 16 cores). 

An example batch script ``my_parallel_job.sh`` is: ::

	#!/bin/bash
	#$ -l rmem=2G
	#$ -pe smp 12
	module load apps/matlab/2017b/binary
	#Run parallel_example.m
	matlab -nodisplay -r parallel_example

where ``parallel_example.m`` is: ::

	%create parallel pool of workers on the local node
	%Ensure that this is the same number as what you requested from the scheduler
	pool = parpool('local',12)
	disp('serial time')
	tic
	n = 200;
	A = 500;
	a = zeros(n);
	for i = 1:n
		a(i) = max(abs(eig(rand(A))));
	end
	toc

	disp('Parallel time')
	tic
	n = 200;
	A = 500;
	a = zeros(n);
	parfor i = 1:n
		a(i) = max(abs(eig(rand(A))));
	end
	toc

	delete(pool)

Parallel MATLAB - Multiple-nodes
--------------------------------

Parallel Matlab using multiple nodes is restricted to 32 cores. 

The user must configure Matlab first by running Matlab interactively and configuring for cluster usage.

This is done by logging into ShARC, launching a qrshx session, module load apps/matlab/2017b & launching matlab. The following command is typed into the command line in the GUI: ::

	configCluster;

Matlab GUI can now be closed.

An example batch script ``submit_Matlab_mpi.sh`` is: ::

	#!/bin/bash
	#$ -M user@sheffield.ac.uk
	#$ -m bea
	#$ -V
	#$ -cwd
	module load apps/matlab/2017b/binary
	#Run parallel_example.m
	matlab -nodisplay -nosplash -r submit_matlab_fnc

where ``submit_matlab_fnc.m`` is: ::

	function submit_matlab_fnc

	cd path_working_directory;
	c=parcluster;
	c.AdditionalProperties.EmailAddress = 'user@sheffield.ac.uk';
	%configure runtime e.g. 40 minutes
	c.AdditionalProperties.WallTime = '00:40:00';
	%configure rmem per process e.g. 4 Gb
	c.AdditionalProperties.AdditionalSubmitArgs = ' -l rmem=4G';
	%parallel_example.m contains the parfor loop, no_of_cores < 31
	j=c.batch(@parallel_example,1,{},'Pool',no_of_cores);

where ``parallel_example.m`` is: ::

	disp('serial time')
	tic
	n = 200;
	A = 500;
	a = zeros(n);
	for i = 1:n
		a(i) = max(abs(eig(rand(A))));
	end
	toc

	disp('Parallel time')
	tic
	n = 200;
	A = 500;
	a = zeros(n);
	parfor i = 1:n
		a(i) = max(abs(eig(rand(A))));
	end
	toc

Note that for multi-node parallel Matlab the maximum number of workers allowed is 31 since the master process requires a parallel licence. Task arrays are supported by all versions, however it is recommended that 2017a (or later) is used 

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
#. Log on to Mathworks site to download the MATLAB installer package for 64-bit Linux ( for R2017b this was called ``matlab_R2017b_glnxa64.zip`` )
#. ``unzip`` the installer package in a directory with ~10GB of space (needed as many MATLAB *archive* files will subsequently be downloaded here).  Using a directory on an NFS mount (e.g. ``/data/${USER}/MathWorks/R2017b``) allows the same downloaded archives to be used to install MATLAB on multiple clusters.
#. ``./install`` to start the graphical installer (needed to download the MATLAB archive files).
#. Select install choice of *Log in to Mathworks Account* and log in with a *License Administrator* account (not a *Licensed End User* (personal) account).
#. Select *Download only*.
#. Select the offered default *Download path* and select the directory you ran ``./install`` from.  Wait a while for all requested archive files to be downloaded.
#. Next, ensure ``installer_input.txt`` looks like the following ::
    
    fileInstallationKey=XXXXX-XXXXX-XXXXX-XXXXX-XXXXX-XXXXX-XXXXX-XXXXX-XXXXX-XXXXX-XXXXX-XXXXX-XXXXX-XXXXX-XXXXX-XXXXX-XXXXX-XXXXX-XXXXX-XXXXX-XXXXX
    agreeToLicense=yes
    outputFile=matlab_2017b_install.log
    mode=silent
    licensePath=/usr/local/packages/matlab/network.lic
    lmgrFiles=false
    lmgrService=false

#. Create the installation directory e.g.: ::

    mkdir -m 2755 -p /usr/local/packages/apps/matlab/R2017b/binary
    chown ${USER}:app-admins /usr/local/packages/apps/matlab/R2017b/binary

#. Run the installer using our customized ``installer_input.txt`` like so: ``./install -mode silent -inputFile ${PWD}/installer_input.txt`` ; installation should finish with exit status ``0`` if all has worked.
#. Ensure the contents of the install directory and the modulefile are writable by those in ``app-admins`` group e.g.: ::

    chmod -R g+w ${USER}:app-admins /usr/local/packages/apps/matlab/R2017b /usr/local/modulefiles/apps/matlab/2017b

The modulefile is
:download:`/usr/local/modulefiles/apps/matlab/2017b/binary </sharc/software/modulefiles/apps/matlab/2017b/binary>`.
