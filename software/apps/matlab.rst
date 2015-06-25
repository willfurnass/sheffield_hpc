MATLAB
======

.. sidebar:: MATLAB 
   
   :Versions:  2013a , 2013b , 2014a, 2015a
   :Support Level: FULL 
   :Dependancies: None
   :URL: http://uk.mathworks.com/products/matlab 
   :Local URL:  http://www.shef.ac.uk/wrgrid/software/matlab
   :Documentation: http://uk.mathworks.com/help/matlab

Scientific Computing and Visualisation 

Interactive Usage
-----------------
After connecting to iceberg (see :ref:`ssh`),  start an interactive sesssion with the :code:`qsh` command.

The lastest version of MATLAB (currently 2015a) is made available with the command

.. code-block:: none

        module load apps/matlab

Alternatively, you can load a specific version with one of of the following commands

.. code-block:: none

       module apps/matlab/2013a
       module apps/matlab/2013b
       module apps/matlab/2014a

You can then run MATLAB by entering :code:`matlab` 

Serial (one CPU) Batch usage
----------------------------
Here, we assume that you wish to run the program :code:`hello.m` on the system.

First, you need to write a batch submission file. We assume you'll call this :code:`my_job.sge` ::

    #!/bin/bash
    #$ -l mem=4G                       # Requiest 4 Gigabytes of memory
    $ -cwd                             # Run job from current directory
    module load apps/matlab            # Make latest version of MATLAB available

    matlab -nosplash -nodisplay -nodesktop -r 'hello'

Ensuring that :code:`hello.m` and :code:`myjob.sge` are both in your current working directory, submit your job to the batch system ::

    qsub my_job.sge

Some notes about this example:

* We are running the script :code:`hello.m` but we drop the `.m` in the call to MATLAB. That is, we do :code:`-r 'hello'` rather than :code:`-r hello.m`.
* All of the :code:`module` commands introducted in the Interactive usage section will also work in batch mode. This allows you to select a specific version of MATLAB if you wish.

Installation notes
------------------
These notes are primarily for system administrators.

Requires the floating license server licserv4.shef.ac.uk to serve the licenses 
for the version of MATLAB to be installed ( or higher versions ) .
An install script and associated files are downloadable from Mathworks site along with all the required toolbox specific installation files. 


