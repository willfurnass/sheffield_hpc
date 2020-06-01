Picard
======

.. sidebar:: Picard

   :Version: 1.129
   :URL: https://github.com/broadinstitute/picard/

A set of Java command line tools for manipulating high-throughput sequencing (HTS) data and formats.
Picard is implemented using the HTSJDK Java library,
supporting accessing of common file formats,
such as SAM and VCF,
used for high-throughput sequencing data.

Interactive Usage
-----------------

After connecting to Iceberg (see :ref:`ssh`), :ref:`start an interactive session <sched_interactive>`
You are likely to want to request at least 8 GB of memory for your session.

Next, load a specific version of Picard using one of the following: ::

   module load apps/binapps/picard/1.129
   module load apps/binapps/picard/1.101

These module commands also update your session's environment to provide Java 1.6 since this is required by Picard 1.129.
An environment variable called ``PICARDHOME`` is created by the module command that contains the path to the requested version of Picard.

You can then run Picard using: ::

   java -jar $PICARDHOME/picard.jar

Installation notes
------------------

**Version 1.129**

A binary install was used.
The binary came from the releases page of the project's GitHub repo: ::

   unzip picard-tools-1.129.zip
   mkdir -p /usr/local/packages6/apps/binapps/picard
   mv ./picard-tools-1.129 /usr/local/packages6/apps/binapps/picard/1.129

**Version 1.101**

A binary install was used.
The binary came from the `project's Sourceforge site <https://sourceforge.net/projects/picard/files/picard-tools/>`__. ::

   unzip picard-tools-1.101.zip
   mv ./picard-tools-1.101 /usr/local/packages6/apps/binapps/picard/1.101/

Modulefile
----------

**Version 1.129**

The module file is on the system at ``/usr/local/modulefiles/apps/binapps/picard/1.129``

Its contents are ::

   #%Module1.0#####################################################################
   ##
   ## Picard 1.129 modulefile
   ##
 
   ## Module file logging
   source /usr/local/etc/module_logging.tcl
   ##
 
   #This version of Picard needs Java 1.6
   module load apps/java/1.6
 
   proc ModulesHelp { } {
           puts stderr "Makes Picard 1.129 available"
   }
 
 
   set version 1.129
   set PICARD_DIR /usr/local/packages6/apps/binapps/picard/$version
 
   module-whatis   "Makes Picard 1.129 available"
 
   prepend-path PICARDHOME $PICARD_DIR

**Version 1.101**

The module file is on the system at ``/usr/local/modulefiles/apps/binapps/picard/1.101``

Its contents are ::

   #%Module1.0#####################################################################
   ##
   ## Picard 1.101 modulefile
   ##
 
   ## Module file logging
   source /usr/local/etc/module_logging.tcl
   ##
 
   #This version of Picard needs Java 1.6
   module load apps/java/1.6
 
   proc ModulesHelp { } {
           puts stderr "Makes Picard 1.101 available"
   }
 
 
   set version 1.101
   set PICARD_DIR /usr/local/packages6/apps/binapps/picard/$version
 
   module-whatis   "Makes Picard 1.101 available"
 
   prepend-path PICARDHOME $PICARD_DIR
