.. _lammps_sharc:

LAMMPS
======

.. sidebar:: LAMMPS

   :Versions:  22_08_2018
   :Support Level: 
   :Dependancies: gcc/4.9.4, openmpi/2.0.1
   :URL: https://lammps.sandia.gov/
   :Documentation: https://lammps.sandia.gov/doc/Manual.html

LAMMPS is a classical molecular dynamics code with a focus on materials modelling. It's an acronym for Large-scale Atomic/Molecular Massively Parallel Simulator.

LAMMPS has potentials for solid-state materials (metals, semiconductors) and soft matter (biomolecules, polymers) and coarse-grained or mesoscopic systems. It can be used to model atoms or, more generically, as a parallel particle simulator at the atomic, meso, or continuum scale.

LAMMPS runs on single processors or in parallel using message-passing techniques and a spatial-decomposition of the simulation domain. Many of its models have versions that provide accelerated performance on CPUs, GPUs, and Intel Xeon Phis. The code is designed to be easy to modify or extend with new functionality.

Interactive Usage
-----------------
After connecting to ShARC (see :ref:`ssh`),  start an interactive session with the ``qrshx`` command.

The latest version of LAMMPS (currently 22_08_2018) is made available by running: ::

   .. code-block:: bash

      module load apps/lammps

Alternatively, you can load a specific version with the following command:

   .. code-block:: bash

      module load apps/lammps/22_08_2018/gcc-4.9.4

You can then run LAMMPS by entering ``lmp``

   .. code-block:: bash

      cp /usr/local/packages/apps/lammps/22_08_2018/gcc-4.9.4/examples.tar.gz .
      tar -xvzf examples.tar.gz
      cd examples/indent # select indent example
      lmp -in in.indent # run indent example
      

Serial (one core) Batch usage
-----------------------------
Here, we assume that you wish to run the program ``helloworld.m`` on the system: ::
	
	function helloworld
		disp('Hello World!')
	end	

First, you need to write a batch submission file. We assume you'll call this ``my_job.sge``: ::

	#!/bin/bash
	#$ -l rmem=4G                  		# Request  4 GB of real memory
	#$ -cwd                        		# Run job from current directory
	module load apps/matlab/2018b/binary  	# Make specific version of MATLAB available
  
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

	# Load the matlab 2018b module
	module load apps/matlab/2018b/binary  

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
	module load apps/matlab/2018b/binary
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

This is done by logging into ShARC, launching a qrshx session, module load apps/matlab/2018a & launching matlab. The following command is typed into the command line in the GUI: ::

	configCluster;

Matlab GUI can now be closed.

An example batch script ``submit_Matlab_mpi.sh`` is: ::

	#!/bin/bash
	#$ -M user@sheffield.ac.uk
	#$ -m bea
	#$ -V
	#$ -j y
	module load apps/matlab/2018b/binary
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
	
	function time = parallel_example
	cd path_working_directory;
	outfile = ['output.txt'];
	fileID = fopen(outfile,'w');
	%disp('Parallel time')
	tic
	n = 200;
	A = 500;
	a = zeros(n);
	parfor i = 1:n
		a(i) = max(abs(eig(rand(A))));
	end
	time=toc;
	fprintf(fileID, '%d', time);
	fclose(fileID);

Note that for multi-node parallel Matlab the maximum number of workers allowed is 31 since the master process requires a parallel licence. Task arrays are supported by all versions, however it is recommended that 2017a (or later) is used. 

MATLAB Engine for Python
------------------------

This is a MathWorks-developed way of running MATLAB from Python.
On ShARC the recommended way of installing this is into a :ref:`conda environment <sharc-python-conda>`.
Here's how you can install the R2017b version into a new conda environment called ``my-environment-name``: ::

    module load apps/python/conda
    conda create -n my-environment-name python=2.7
    source activate my-environment-name 

    pushd /usr/local/packages/apps/matlab/2017b/binary/extern/engines/python
    python setup.py build -b $TMPDIR install
    popd

`More information <https://uk.mathworks.com/help/matlab/matlab_external/install-the-matlab-engine-for-python.html>`__ on the MATLAB Engine for Python,
including basic usage.

Training
--------

* CiCS run an `Introduction to Matlab course <http://rcg.group.shef.ac.uk/courses/matlab/>`_
* In November 2015, CiCS hosted a masterclass in *Parallel Computing in MATLAB*. The materials `are available online <http://rcg.group.shef.ac.uk/courses/mathworks-parallelmatlab/>`_


Installation notes
------------------

These notes are primarily for system administrators.

Installation and configuration is a five-stage process:

* Set up the floating license server (the license server for earlier MATLAB versions can be used), ensuring that it can serve licenses for any new versions of MATLAB that you want to install
* Run a graphical installer to download MATLAB *archive* files used by the main (automated) installation process
* Run the same installer in 'silent' command-line mode to perform the installation using those archive files and a text config file.
* Install a relevant modulefile
* Configure MATLAB parallel (multi-node)

In more detail:

#. If necessary, update the floating license keys on ``licserv4.shef.ac.uk`` to ensure that the licenses are served for the versions to install.
#. Log on to Mathworks site to download the MATLAB installer package for 64-bit Linux ( for R2018b this was called ``matlab_R2018b_glnxa64.zip`` )
#. ``unzip`` the installer package in a directory with ~10GB of space (needed as many MATLAB *archive* files will subsequently be downloaded here).  Using a directory on an NFS mount (e.g. ``/data/${USER}/MathWorks/R2018a``) allows the same downloaded archives to be used to install MATLAB on multiple clusters.
#. ``./install`` to start the graphical installer (needed to download the MATLAB archive files).
#. Select install choice of *Log in to Mathworks Account* and log in with a *License Administrator* account (not a *Licensed End User* (personal) account).
#. Select *Download only*.
#. Select the offered default *Download path* and select the directory you ran ``./install`` from.  Wait a while for all requested archive files to be downloaded.
#. Next, ensure ``installer_input.txt`` looks like the following ::
    
    fileInstallationKey=XXXXX-XXXXX-XXXXX-XXXXX-XXXXX-XXXXX-XXXXX-XXXXX-XXXXX-XXXXX-XXXXX-XXXXX-XXXXX-XXXXX-XXXXX-XXXXX-XXXXX-XXXXX-XXXXX-XXXXX-XXXXX
    agreeToLicense=yes
    outputFile=matlab_2018b_install.log
    mode=silent
    licensePath=/usr/local/packages/matlab/network.lic
    lmgrFiles=false
    lmgrService=false

#. Create the installation directory e.g.: ::

    mkdir -m 2755 -p /usr/local/packages/apps/matlab/R2018b/binary
    chown ${USER}:app-admins /usr/local/packages/apps/matlab/R2018b/binary

#. Run the installer using our customized ``installer_input.txt`` like so: ``./install -mode silent -inputFile ${PWD}/installer_input.txt`` ; installation should finish with exit status ``0`` if all has worked.
#. Ensure the contents of the install directory and the modulefile are writable by those in ``app-admins`` group e.g.: ::

    chmod -R g+w ${USER}:app-admins /usr/local/packages/apps/matlab/R2018b /usr/local/modulefiles/apps/matlab/2018b

#. The modulefile is :download:`/usr/local/modulefiles/apps/matlab/2018b/binary </sharc/software/modulefiles/apps/matlab/2018b/binary>`.

#. Copy integration scripts to MATLAB local directory (required for MATLAB parallel (multi-node)): ::

    cd /usr/local/packages/apps/matlab/2018b/binary/toolbox/local
    cp -r /usr/local/packages/apps/matlab/parallel_mpi_integration_scripts/* .

#. R2018a Update 4 to mitigate Matlab crashes on Centos 7.5. Download R2018a Update 4 installer. Copy to ShARC, and run using ./R2018a_Update_4_glnxa64.sh, and specify install directory /usr/local/packages/matlab/2018a/binary
