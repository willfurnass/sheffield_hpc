Annovar
=======

.. sidebar:: Annovar

   :Version:  2015DEC14
   :URL: http://annovar.openbioinformatics.org/en/latest/

Cufflinks assembles transcripts, estimates their abundances, and tests for differential expression and regulation in RNA-Seq samples.
It accepts aligned RNA-Seq reads and assembles the alignments into a parsimonious set of transcripts. Cufflinks then estimates the relative abundances of these transcripts based on how many reads support each one, taking into account biases in library preparation protocols.


Interactive Usage
-----------------
After connecting to Iceberg (see :ref:`ssh`),  start an interactive session with the :code:`qsh` or `qrsh` command.

To add the annovar binaries to the system PATH, execute the following command ::

        module load apps/binapps/annovar/2015DEC14

Documentation
-------------
The annovar manual is available online at http://annovar.openbioinformatics.org/en/latest/

Installation notes
------------------
The install is a collection of executable Perl scripts. Installation involves copying the directory and adding it to the PATH ::
It seems that annovar uses dates to distinguish between releases rather than version numbers

    tar -xvzf ./annovar.latest.tar.gz
    mkdir -p mkdir -p /usr/local/packages6/apps/binapps/annovar
    mv annovar /usr/local/packages6/apps/binapps/annovar/2015DEC14

Modulefile
----------
* The module file is on the system at `/usr/local/modulefiles/apps/binapps/annovar/2015DEC14`

Its contents are ::

  #%Module1.0#####################################################################
  ##
  ## annovar modulefile
  ##

  ## Module file logging
  source /usr/local/etc/module_logging.tcl
  ##

  proc ModulesHelp { } {
          puts stderr "Makes annovar 2015DEC14 available"
  }

  set ANNOVAR_DIR /usr/local/packages6/apps/binapps/annovar/2015Dec14

  module-whatis   "Makes annovar 2015DEC14 available"

  prepend-path PATH $ANNOVAR_DIR
