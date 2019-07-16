.. _matlab_sharc:

MATLAB
======

.. sidebar:: MATLAB

   :Versions:  2019a
   :Support Level: FULL
   :Dependancies: None
   :URL: http://uk.mathworks.com/products/matlab
   :Documentation: http://uk.mathworks.com/help/matlab

Scientific computing and visualisation.


Interactive usage
-----------------
After connecting to Bessemer,  start an interactive session with the ``srun --pty bash –i`` command.

The latest version of MATLAB (currently 2019a) is made available by running:

.. code-block:: bash

   module use /usr/local/modulefiles/staging/eb/all/
   module load MATLAB/2019a/binary

You can then run MATLAB by entering ``matlab &``.


Serial (one CPU) batch usage
----------------------------
Here, we assume that you wish to run the program ``helloworld.m`` on the system:
	
.. code-block:: matlab

   function helloworld
       disp('Hello World!')
   end	

First, you need to write a batch submission file.
We assume you'll call this ``my_job.sh``:

.. code-block:: bash

   #!/bin/bash
   #SBATCH --mem=16000                # Request  16 GB of real memory

   module use /usr/local/modulefiles/staging/eb/all/
   module load MATLAB/2019a/binary

   matlab -nodesktop -nosplash -r helloworld

Ensure that ``helloworld.m`` and ``my_job.sh`` are both in your current working directory, 
then submit your job to the batch system:

.. code-block:: bash

   sbatch my_job.sge

Note that we are running the script ``helloworld.m`` 
but we drop the ``.m`` in the call to MATLAB. 
That is, we do ``-r helloworld`` 
rather than ``-r helloworld.m``. 
The output will be written to the job text file when the job finishes.


Parallel MATLAB
---------------

Parallel MATLAB using multiple nodes is restricted to 40 cores. 

The user must first configure MATLAB for cluster usage by starting MATLAB interactively.
This is done by logging into ShARC, 
launching a ``qrshx`` session, 
loading a version of MATLAB (e.g. using ``module load apps/matlab/2019a``) and 
launching MATLAB with ``matlab``. 
You then need to type the following at the prompt within the MATLAB GUI:


An example (using 40 cores) batch script ``submit_Matlab_par.sh`` is:

.. code-block:: bash

   #!/bin/bash
   #SBATCH --nodes=1
   #SBATCH --mem=16000
   #SBATCH --ntasks-per-node=40
   #SBATCH --time=00:05:00
   #SBATCH --job-name=matlab_batch_test

   module use /usr/local/modulefiles/staging/eb/all/
   module load MATLAB/2019a/binary

   matlab -nodisplay -nosplash -r parallel_example


where ``parallel_example.m`` is:

.. code-block:: matlab

   outfile = ['output.txt'];
   fileID = fopen(outfile,'w');
   pool = parpool('local',16)
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


Note: parallel_example.m creates 200 square (500x500) matrices comprised of random values and calculates the eigenvalues of each (and records the maximum eigenvalue for each matrix in the array a).
