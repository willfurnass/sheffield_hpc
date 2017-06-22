Abaqus
======

.. sidebar:: Abaqus
   
   :Version: 6.14-2
   :Dependencies: Intel compiler 15.0.7 and Foxit modules loaded. User subroutines need Intel compiler 2011 or above, GCC 4.1.2 or above. 
   :URL: http://www.3ds.com/products-services/simulia/products/abaqus/ 
   :Documentation: https://www.shef.ac.uk/wrgrid/software/abaqus


Abaqus is a software suite for Finite Element Analysis (FEA) developed by Dassault Systèmes.


Usage
-----

Abaqus 6.14-2 can be activated using the module file::

    module load apps/abaqus/614
	
Type ``abaqus cae`` to launch the Abaqus GUI during an interactive session with X Window support (e.g. an interactive ``qsh`` session).
Type ``abaqus`` for the command line interface. Typing ``abaqus -help`` will display a list of usage options.
The PDF viewer ``foxit`` can be launched to view the PDF documentation located at "/usr/local/packages/apps/abaqus/Documentation/docs/v6.14/pdf_books".


Abaqus example problems
-----------------------

Abaqus contains a large number of example problems which can be used to become familiar with Abaqus on the system.
These example problems are described in the Abaqus documentation and can be obtained using the Abaqus ``fetch`` command.
For example, after loading the Abaqus module enter the following at the command line to extract the input file for test problem s4d::

    abaqus fetch job=s4d
	
This will extract the input file ``s4d.inp`` to run the computation defined by the commands and batch submission script below.


Batch jobs
----------

The easiest way of running a batch job for a particular version of Abaqus (e.g. 6.14-2) is::
    
    module load apps/abaqus/614
    runabaqus
	
The ``runabaqus`` command submits an Abaqus input file into the batch queuing system and can take a number of different parameters according to your requirements.
Typing ``runabaqus`` will display information about its usage.

Users are encouraged to write their own batch submission scripts. The following is an example batch submission script, ``my_job.sh``, to run the executable ``abq6142`` and which is submitted to the queue by typing ``qsub my_job.sh``::

    #!/bin/bash
    #$ -cwd
    #$ -l h_rt=00:30:00
    #$ -l rmem=2G
    #$ -pe smp 4

    module load apps/abaqus/614

    mkdir /scratch/$USER

    abq6142 job=my_job input=s4d.inp scratch=/scratch/$USER memory="2gb" interactive mp_mode=threads cpus=$NSLOTS
	
The above script requests 4 cores using the OpenMP parallel environment ``smp`` with a runtime of 30 mins and 2G of real memory per core. The Abaqus input file is ``s4d.inp``.

*User subroutines:* The script below is an example of a batch submission script for a single core job with a runtime of 30 mins, 8G of real memory and with user subroutine ``umatmst3.f`` and input file ``umatmst3.inp``::

    #!/bin/bash
    #$ -cwd
    #$ -l h_rt=00:30:00
    #$ -l rmem=8G

    module load apps/abaqus/614
    
    mkdir /scratch/$USER
    
    abq6142 job=my_user_job input=umatmst3.inp user=umatmst3.f scratch=/scratch/$USER memory="8gb" interactive

The input file ``umatmst3.inp`` and the Fortran user subroutine ``umatmst3.f`` are obtained by typing ``abaqus fetch job=umatmst3*``.
Note that the module ``dev/intel-compilers/15.0.7``, required for user subroutines, is automatically loaded when the module for Abaqus is loaded.  

**Important information:** Please note that at present Abaqus will not run on more than one node when using MPI on ShARC.


Installation notes
------------------

Abaqus 6.14-2 was installed using the
:download:`install_abaqus.sh </sharc/software/install_scripts/apps/abaqus/614/install_abaqus.sh>` script; the module
file is
:download:`614 </sharc/software/modulefiles/apps/abaqus/614>`. The binary installation was tested by launching ``abaqus cae`` and by using the above batch submission scripts.
Abaqus at present does not run on more than one node when using MPI due to password-less ssh being disabled across nodes on ShARC.
