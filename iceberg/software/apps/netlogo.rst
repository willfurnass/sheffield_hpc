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

Installation notes
------------------
Download and untar the installer to `/usr/local/packages6/apps/binapps/netlogo/netlogo-5.3.1-64`

Modulefile
----------
* The module file is on the system at `/usr/local/modulefiles/apps/binapps/netlogo/5.3.1`
