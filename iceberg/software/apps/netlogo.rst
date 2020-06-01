.. _netlogo:

Netlogo
=======

.. sidebar:: Netlogo

   :Version:  5.3.1
   :URL: https://ccl.northwestern.edu/netlogo/index.shtml

NetLogo is a multi-agent programmable modeling environment.
It is used by tens of thousands of students, teachers and researchers worldwide.

Interactive Usage
-----------------
After connecting to Iceberg (see :ref:`ssh`),  :ref:`start an interactive session <sched_interactive>` then
load a specific version of NetLogo using: ::

   module load apps/binapps/netlogo/5.3.1

This command makes the NetLogo executables available to your session by adding the install directory to your ``PATH`` variable.
Start the NetLogo Graphical User Interface with the command: ::

   NetLogo

Batch usage
-----------

Here is an example batch submission file for a Netlogo model that contains a `BehaviourSpace <https://ccl.northwestern.edu/netlogo/docs/behaviorspace.html>`_ experiment called **ParallelExperiment1**.
Using BehaviourSpace is the easiest way to exploit parallel processors in Netlogo. ::


   #!/bin/bash
   # Request 3 Gigabytes of RAM per core
   #$ -l rmem=3G
   #Combine stderr and stdout into one file
   #$ -j y
   #Make sure that the value below matches the number of threads
   #$ -pe openmp 4
   #Make sure this matches the number of openmp slots requested
   threads=4
 
   module load apps/binapps/netlogo/5.3.1
 
   #Set the filename of the output csv file
   #Will be called out_4.csv in this case
   output_table=out_$threads.csv
   #Experiment and model names
   experiment=ParallelExperiment1
   model=ParallelPredation.nlogo
 
   echo "$threads threads requested"
 
   #You have to run netlogo from its install directory in order for extensions to work
   cd $NETLOGO_DIR
   java -Xmx1024m -Dfile.encoding=UTF-8 -cp $NETLOGO_DIR/NetLogo.jar org.nlogo.headless.Main --model $SGE_O_WORKDIR/$model --experiment $experiment --table $SGE_O_WORKDIR/$output_table --threads $threads

Call the above file ``submit_netlogo.sh`` and
submit it with ``qsub submit_netlogo.sh``.

All files required to run this example can be found in the `HPC Examples <https://github.com/mikecroucher/HPC_Examples>`_ GitHub repository.

Netlogo uses shared memory parallelism.
As such, the maximum number of OpenMP slots you can request on Iceberg is 16 (the maximum number of CPU cores per node).

Troubleshooting
---------------
**Problems with X Windows**

When you run the NetLogo command, you get the following error: ::

   NetLogo Error invoking method.
   NetLogo Failed to launch JVM

This could be because you have have :ref:`started a non-graphical interactive session rather than a graphical one <sched_interactive>` and/or
have not enabled X forwarding if you connected using :ref:`SSH <ssh>`.

Installation notes
------------------
Download and untar the installer to ``/usr/local/packages6/apps/binapps/netlogo/netlogo-5.3.1-64``.

**Post-installation actions**

One NetLogo user required the Dynamic scheduler extension.
Installed by doing ::

   wget https://github.com/colinsheppard/Dynamic-Scheduler-Extension/archive/v0.2.1.tar.gz
   tar -xvzf ./v0.2.1.tar.gz
   mv Dynamic-Scheduler-Extension ./dynamic-scheduler
   mv ./dynamic-scheduler/ /usr/local/packages6/apps/binapps/netlogo/netlogo-5.3.1-64/app/extensions/

Modulefile
----------

The module file is on the system at ``/usr/local/modulefiles/apps/binapps/netlogo/5.3.1``
