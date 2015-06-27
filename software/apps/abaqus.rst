Abaqus
======

.. sidebar:: Abaqus

   :Versions:  6.13,6.12 and 6.11
   :Support Level: FULL
   :Dependancies: None
   :URL: http://www.3ds.com/products-services/simulia/products/abaqus/
   :Local URL:  https://www.shef.ac.uk/wrgrid/software/abaqus

Abaqus is a software suite for Finite Element Analysis (FEA) developed by Dassault Systèmes.

Interactive usage
-----------------
After connecting to iceberg (see :ref:`ssh`),  start an interactive sesssion with the :code:`qsh` command. Alternatively, if you require more memory, for example 16 gigabytes, use the command :code:`qsh -l mem=16G` 

The lastest version of Abaqus (currently version 6.13) is made available with the command

.. code-block:: none

        module load apps/abaqus

Alternatively, you can make a specific version available with one of the following commands

.. code-block:: none

      module load apps/abaqus/613
      module load apps/abaqus/612 
      module load apps/abaqus/611 

After that, simply type :code:`abaqus` to get the command-line interface to abaqus or type :code:`abaqus cae` to get the GUI interface.

Abaqus example problems
-----------------------
Abaqus contains a large number of example problems which can be used to become familiar with Abaqus on the CSF. These example problems are described in the Abaqus documentation, and can be obtained using the Abaqus fetch command. For example, after loading the Abaqus module enter the following at the command line to extract the input file for test problem s4d ::

    abaqus fetch job=s4d

This will extract the input file s4d.inp, to run the computation defined by this input file replace input=myabaqusjob with input=s4d in the commands and scripts below.

Batch submission of a single core job
-------------------------------------
In this example, we will run the s4d.inp file on a single core using 8 Gigabytes of memory.  After connecting to iceberg (see :ref:`ssh`),  start an interactive sesssion with the :code:`qrsh` command. 

Load the latest version of Abaqus and fetch the s4d example ::

    module load apps/abaqus
    abaqus fetch job=s4d

Now, you need to write a batch submission file. We assume you'll call this :code:`my_job.sge` ::

    #!/bin/bash
    #$ -S /bin/bash
    #$ -cwd
    #$ -l rmem=8G
    #$ -l mem=8G
    
    module load apps/abaqus    

    abaqus job=my_job input=s4d.inp scratch=/scratch memory="8gb" interactive

Submit the job with the command ``qsub my_job.sge``

Important notes:

* We have requested 8 gigabytes of memory in the above job. The ``memory="8gb"`` switch tells abaqus to use 8 gigabytes. The ``#$ -l rmem=8G`` and ``#$ -l mem=8G`` tells the system to reserve 8 gigabytes of real and virtual memory resptively. It is important that these numbers match.
* Note the word ``interactive`` at the end of the abaqus command. Your job will not run without it.



