.. _ansys:

Ansys
=====

.. sidebar:: Ansys

   :Version:  14 , 14.5 , 15
   :Support Level: FULL
   :Dependencies: If using the User Defined Functions (UDF) will also need the following:
                  For Ansys Mechanical, Workbench, CFX and AutoDYN : INTEL 14.0 or above Compiler
                  For Fluent :  GCC 4.6.1 or above 
   :URL: http://www.ansys.com/en_uk
   :Local URL: http://www.shef.ac.uk/cics/research/software/fluent

ANSYS suite of programs can be used to numerically simulate a large variety of structural and fluid dynamics problems found in many engineering, physics, medical, aeronotics and automative industry applications

Interactive usage
-----------------
After connecting to iceberg (see :ref:`ssh`),  start an interactive sesssion with the :code:`qsh` command. Alternatively, if you require more memory, for example 16 gigabytes, use the command :code:`qsh -l rmem=16G` 

To make the latest version of Ansys available, run the following module command

.. code-block:: none

        module load apps/ansys

Alternatively, you can make a specific version available with one of the following commands

.. code-block:: none

      module load apps/ansys/14
      module load apps/ansys/14.5 
      module load apps/ansys/15.0 

After loading the modules, you can issue the following commands to run an ansys product. The teaching versions mentioned below will only allow models for upto 500,000 elements but because we have plenty of teaching licenses it will always be possible to use without the fear of running out of licenses. 

.. code-block:: none

      ansyswb  : to run Ansys workbench      
      ansys    : to run Ansys Mechanical outside the workbench
      ansystext: to run line-mode version of ansys
      Ansys and Ansystext : to run the teaching license version of the above two commands.
      fluent   : to run Fluent outside the workbench
      fluentext: to run the line-mode version of Fluent
      Fluent and Fluentext : to run the teaching license version of the above two commands.
      icemcfx or icem: to run icemcfd outside the workbench.
       


Running Batch fluent and ansys jobs
-----------------------------------


The easiest way of running batch ansys and fluent jobs is to select and use one of the the following commands according to the version of fluent or ansys you wish to run:

.. code-block:: none

     runfluent14 , runfluent145 , runfluent150 ,  runansys14 , runansys145 , runansys150 




