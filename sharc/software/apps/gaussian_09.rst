Gaussian 09
===========

.. sidebar:: Gaussian 09

   :Version: D.01
   :Dependencies: Portland Group Fortran compiler and GCC compiler. Module for PGI 17.5 loaded.
   :URL: http://gaussian.com
   :Documentation: http://gaussian.com/techsupport/


Gaussian 09 provides state-of-the-art capabilities for electronic structure modelling. GaussView is the graphical interface used with Gaussian.


Usage
-----

Gaussian 09 Revision D.01 and GaussView 5.0.8, the GUI for Gaussian 09, can be activated using the module file::

    module load apps/gaussian_09/d.01/pgi-17.5

The GaussView executable is ``gv``. ``gv`` can be launched from an interactive session with X Window support (e.g. an interactive ``qsh`` session).
The Gaussian 09 command line executuable is ``g09``.

**Important note:** Access to Gaussian 09 on ShARC is available to all users from the Department of Chemistry by default (i.e. user IDs beginning with ``ch``);
other users wishing to access Gaussian 09 will need to contact ``research-it@sheffield.ac.uk`` to request access (i.e. to be added to unix group ``gaussian09``).


Gaussian 09 utilities
---------------------

The Gaussian 09 utilities can be accessed by loading the Gaussian 09 module file (see above).
The utilities are executed by typing their name, e.g., ``formchk`` and ``cubegen``. The full list of utilities is described at http://gaussian.com/utils.


Batch jobs
----------

Users are encouraged to write their own batch submission scripts. The following is an example batch submission script, ``my_job.sh``, to run ``g09`` and which is submitted to the queue by typing ``qsub my_job.sh``. ::

    #!/bin/bash
    #$ -cwd
    #$ -l h_rt=00:30:00
    #$ -l rmem=2G
    #$ -pe smp 2

    module load apps/gaussian_09/d.01/pgi-17.5
    
    export GAUSS_SCRDIR=$TMPDIR

    g09 my_input.gjf

The script requests 2 cores using the shared memory parallel environment ``smp`` with a runtime of 30 mins and 2GB of real memory per core. The Gaussian 09 input file is ``my_input.gjf``.

**Please note:** The above script specifies that the scratch directory for Gaussian 09, ``GAUSS_SCRDIR``, is set to ``$TMPDIR`` on a compute node.
 

Installation notes
------------------

Gaussian 09 Revision D.01 and GaussView 5.0.8 were installed using the
:download:`install_gaussian_09.sh </sharc/software/install_scripts/apps/gaussian_09/d.01/pgi-17.5/install_gaussian_09.sh>` script;
the module file is
:download:`pgi-17.5 </sharc/software/modulefiles/apps/gaussian_09/d.01/pgi-17.5>`.
The Gaussian 09 code was compiled with PGI 17.5 and GCC 4.8.5. 

The Gaussian 09 installation was tested by running a batch job using the following text (including a blank line at the end) in an input file and the ``my_job.sh`` batch script, above. ::

    %chk=h2o.chk
    %nproc=2
    %mem=2GB
    #n hf/6-31G(d,p) opt freq

    H2O

    0 1
    O
    H 1 r1
    H 1 r2 2 a1

    r1 1.0
    r2 1.0
    a1 105.0

