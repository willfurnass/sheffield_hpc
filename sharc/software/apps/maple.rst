.. _maple_sharc:

Maple
=====

.. sidebar:: Maple

   :Latest Version:  2018
   :URL: http://www.maplesoft.com/products/maple/

Scientific Computing and Visualisation

Interactive Usage
-----------------
After connecting to ShARC (see :ref:`ssh`), :ref:`start an interactive graphical session <sched_interactive>` then
load a specific version of Maple using: ::

   module load apps/maple/2018/binary

You can then run the graphical version of Maple by entering ``xmaple`` or the command line version by entering ``maple``.

Batch usage
-----------

It is not possible to run Maple worksheets in batch mode.
Instead, you must convert your worksheet to a pure text file that contains a set of Maple input commands.
You can do this in Maple by opening your worksheet and clicking on **File->Export As->Maple Input**.
The result will have the file extension ``.mpl``.

An example SGE submission script that makes use of a ``.mpl`` file called, for example, ``mycode.mpl`` is: ::

   #!/bin/bash
   # Request 4 gigabytes of real memory
   #$ -l rmem=4G

   module load apps/maple/2018

   maple < mycode.mpl

For general information on how to submit batch jobs refer to :ref:`sched_batch`.

Tutorials
---------

* `High Performance Computing with Maple <https://rse.shef.ac.uk/blog/hpc-maple-1/>`_: A tutorial from the Sheffield Research Software Engineering group on how to use Maple in a High Performance Computing environment

Installation notes
------------------

These are primarily for administrators of the system.

**Maple 2018**

* Run the installer **Maple2018.1LinuxX64Installer.run** and follow instructions.
* Choose Install folder ``/usr/local/packages/apps/maple/2018/binary``
* Do not configure MATLAB
* Choose a network license. Details on IT Services internal wiki.
* Uncheck 'Enable periodic checking for Maple 2018 updates'
* Check 'Check for updates now'

The module file is :download:`/usr/local/modulefiles/apps/maple/2018/binary </sharc/software/modulefiles/apps/maple/2018/binary>`.
