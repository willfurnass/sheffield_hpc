.. _blast:

BLAST
=====

.. sidebar:: ncbi-blast

   :Version:  2.3.0
   :URL: https://blast.ncbi.nlm.nih.gov/

BLAST+ is a new suite of BLAST tools that utilizes the NCBI C++ Toolkit.
The BLAST+ applications have a number of performance and feature improvements over the legacy BLAST applications.
For details, please see the `BLAST+ user manual <https://www.ncbi.nlm.nih.gov/books/NBK279690/>`_.

Interactive Usage
-----------------

After connecting to Iceberg (see :ref:`ssh`), :ref:`start an interactive session <sched_interactive>` then
load a specific version of BLAST+ with: ::

   module load apps/binapps/ncbi-blast/2.3.0

This command makes the BLAST+ executables available to your session by adding the install directory to your ``PATH`` variable.
It also sets the ``BLASTDB`` database environment variable.

You can now run commands directly. e.g.: ::

   blastn -help

Databases
---------

The following databases have been installed following a user request

* **nr.*tar.gz** Non-redundant protein sequences from GenPept, Swissprot, PIR, PDF, PDB, and NCBI RefSeq
* **nt.*tar.gz** Partially non-redundant nucleotide sequences from all traditional divisions of GenBank, EMBL, and DDBJ excluding GSS,STS, PAT, EST, HTG, and WGS.

A full list of databases is `available on the NCBI FTP site <https://ftp.ncbi.nlm.nih.gov/blast/documents/blastdb.html>`__.

If you need any of these installing, please contact IT Services.

Installation notes
------------------

This was an install from binaries: ::

   wget ftp://ftp.ncbi.nlm.nih.gov/blast/executables/LATEST/ncbi-blast-2.3.0+-x64-linux.tar.gz
   tar -xzf ./ncbi-blast-2.3.0+-x64-linux.tar.gz
   mkdir -p /usr/local/packages6/apps/binapps/ncbi-blast/
   mv ncbi-blast-2.3.0+ /usr/local/packages6/apps/binapps/ncbi-blast/


   # Create database directory
   mkdir -p /usr/local/packages6/apps/binapps/ncbi-blast/ncbi-blast-2.3.0+/db

   # Install the nr database
   cd /usr/local/packages6/apps/binapps/ncbi-blast/ncbi-blast-2.3.0+/db
   for num in `seq 0 48`;
   do
   paddednum=`printf "%02d" $num`
   `wget ftp://ftp.ncbi.nlm.nih.gov/blast/db/nr.$paddednum.tar.gz`
   done

   # Install the nt database
   for num in `seq 0 36`;
   do
   paddednum=`printf "%02d" $num`
   `wget ftp://ftp.ncbi.nlm.nih.gov/blast/db/nt.$paddednum.tar.gz`
   done

   for f in *.tar.gz; do tar -xvzf $f; done

Testing
-------
No testing has been performed.
If you can suggest a suitable test suite, please contact us.

Modulefile
----------

The module file is on the system at ``/usr/local/modulefiles/apps/binapps/ncbi-blast/2.3.0``.

The contents of the module file is ::

  #%Module1.0#####################################################################
  ##
  ## BLAST 2.3.0 modulefile
  ##

  ## Module file logging
  source /usr/local/etc/module_logging.tcl
  ##

  proc ModulesHelp { } {
          puts stderr "Makes BLAST 2.3.0 available"
  }

  set BLAST_DIR /usr/local/packages6/apps/binapps/ncbi-blast/ncbi-blast-2.3.0+

  module-whatis   "Makes BLAST 2.3.0 available"

  prepend-path PATH $BLAST_DIR/bin
  prepend-path BLASTDB $BLAST_DIR/db
