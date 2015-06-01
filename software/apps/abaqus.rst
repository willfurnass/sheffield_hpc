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

      module load apps/abaqus/611 
      module load apps/abaqus/612 
      module load apps/abaqus/613

After that, simply type :code:`abaqus` to get the command-line interface to abaqus or type :code:`abaqus cae` to get the GUI interface.

