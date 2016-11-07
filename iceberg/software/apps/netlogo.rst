.. _netlogo:

Netlogo
=======

.. sidebar:: Netlogo

   :Version:  5.3.1
   :URL: https://ccl.northwestern.edu/netlogo/index.shtml

NetLogo is a multi-agent programmable modeling environment. It is used by tens of thousands of students, teachers and researchers worldwide.

Interactive Usage
-----------------
After connecting to Iceberg (see :ref:`ssh`),  start an interactive session with the :ref:`qrshx` command.

The latest version of NetLogo, currently 5.3.1, is made available with the command ::

    module load apps/binapps/netlogo/

Alternatively, you can load a specific version with ::

    module load apps/binapps/netlogo/5.3.1

This command makes the NetLogo executables available to your session by adding the install directory to your PATH variable.
Start the NetLogo Graphical User Interface with the command ::

    NetLogo

Troubleshooting
---------------
**Problems with X Windows**

When you run the NetLogo command, you get the following error ::

  NetLogo Error invoking method.
  NetLogo Failed to launch JVM

This could be because you have not requested an X session properly. Ensure that you have logged into the system using the `-X` flag e.g. ::

  ssh -X abcd1@iceberg.sheffield.ac.uk

Also ensure that you have requested a :ref:`qrshx` session rather than a :ref:`qrsh` session.

**Java memory issues**

You get the following error message ::

  Error occurred during initialization of VM
  Could not reserve enough space for object heap
  Error: Could not create the Java Virtual Machine.
  Error: A fatal exception has occurred. Program will exit.

This is the 'Virtual Memory' problem described in the :ref:`Java-iceberg` section. You need to request more memory for your session.
For example, request an interactive session like this ::

    qrshx -l mem=8G -l rmem=8G

Installation notes
------------------
Download and untar the installer to `/usr/local/packages6/apps/binapps/netlogo/netlogo-5.3.1-64`

Modulefile
----------
* The module file is on the system at `/usr/local/modulefiles/apps/binapps/netlogo/5.3.1`
