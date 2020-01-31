Abaqus
======

.. sidebar:: Abaqus
   
   :Versions: 2018, 2017-ga, 2017, 6.14-2
   :Dependencies: Module loaded for Intel compiler 15.0.7 (and Foxit for Abaqus version 6.14-2). User subroutines need Intel compiler 2011 or above, GCC 4.1.2 or above. 
   :URL: http://www.3ds.com/products-services/simulia/products/abaqus/ 
   :Documentation: https://www.3ds.com/support/documentation/users-guides/
   :Local URL: https://www.sheffield.ac.uk/it-services/research/software/abaqus


Abaqus is a software suite for Finite Element Analysis (FEA) developed by Dassault Systèmes.


Usage
-----

Abaqus versions 2018, 2017-ga, 2017 or 6.14-2 can be activated using the module files::

    module load apps/abaqus/2018/binary
    module load apps/abaqus/2017-ga/binary
    module load apps/abaqus/2017/binary
    module load apps/abaqus/6.14-2/binary
	
Type ``abaqus cae`` to launch the Abaqus GUI from an interactive session with X Window support (e.g. an interactive ``qsh`` session). Please see usage note below for graphics support options.
Type ``abaqus`` for the command line interface. Typing ``abaqus -help`` will display a list of usage options.

Abaqus 2017-ga (module ``apps/abaqus/2017-ga/binary``) has the Tosca component installed and is equivalent to Abaqus 2017 ('ga' is an accronym for 'general availabilty').

The PDF viewer ``foxit`` can be launched to view the PDF documentation for Abaqus 6.14-2 located at ``/usr/local/packages/apps/abaqus/6.14-2/binary/Documentation/docs/v6.14/pdf_books``.

Abaqus 2017 documentation in HTML format is located at ``/usr/local/packages/apps/abaqus/2017/binary/Documentation/DSSIMULIA_Established_homepage_English.htm``.

Abaqus 2017-ga documentation in HTML format is located at ``/usr/local/packages/apps/abaqus/2017-ga/binary/SIMULIA2017doc/DSSIMULIA_Established_homepage_English.htm``.

Abaqus 2018 documentation in HTML format is located at ``/usr/local/packages/apps/abaqus/2018/binary/SIMULIA2018doc/DSSIMULIA_Established_homepage_English.htm``.


**Note:** When using hardware-accelerated graphics rendering for Abaqus 6.14-2 on ShARC, i.e., during a ``qsh-vis`` interactive session, please run ``abq6142 cae`` to launch the GUI. When using a general compute node for Abaqus 2017, 2017-ga or 2018 on ShARC, please run ``abaqus cae -mesa`` or ``abq2017 cae -mesa`` (or ``abq2018 cae -mesa``) to launch the GUI without support for hardware-accelerated graphics rendering. The option ``-mesa`` disables hardware-accelerated graphics rendering within Abaqus's GUI.


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
    
    module load apps/abaqus/6.14-2/binary
    runabaqus
	
The ``runabaqus`` command submits an Abaqus input file into the batch queuing system and can take a number of different parameters according to your requirements.
Typing ``runabaqus`` will display information about its usage. **Note:** ``runabaqus`` is not setup for use with Abaqus 2017, 2017-ga and 2018.

Users are encouraged to write their own batch submission scripts. The following is an example batch submission script, ``my_job.sh``, to run the executable ``abq6142`` and which is submitted to the queue by typing ``qsub my_job.sh``::

    #!/bin/bash
    #$ -cwd
    #$ -l h_rt=00:30:00
    #$ -l rmem=2G
    #$ -pe smp 4

    module load apps/abaqus/6.14-2/binary

    abq6142 job=my_job input=s4d.inp scratch=$TMPDIR memory="2gb" interactive mp_mode=threads cpus=$NSLOTS
	
The above script requests 4 cores using the OpenMP parallel environment ``smp`` with a runtime of 30 mins and 2 GB of real memory per core. The Abaqus input file is ``s4d.inp``.

**User subroutines:** The script below is an example of a batch submission script for a single core job with a runtime of 30 mins, 8 GB of real memory and with user subroutine ``umatmst3.f`` and input file ``umatmst3.inp``::

    #!/bin/bash
    #$ -cwd
    #$ -l h_rt=00:30:00
    #$ -l rmem=8G

    module load apps/abaqus/6.14-2/binary
    
    abq6142 job=my_user_job input=umatmst3.inp user=umatmst3.f scratch=$TMPDIR memory="8gb" interactive

The input file ``umatmst3.inp`` and the Fortran user subroutine ``umatmst3.f`` are obtained by typing ``abaqus fetch job=umatmst3*``.
Note that the module ``dev/intel-compilers/15.0.7``, required for user subroutines, is automatically loaded when the module for Abaqus is loaded.  

**Important information:** Please note that at present Abaqus will not run on more than one node when using MPI on ShARC. The SGE option ``-l excl=true`` can be used to request that an MPI job runs on one compute node only. The recommended way to run Abaqus in parallel on ShARC is to use OpenMP.

Using /fastdata as your Abaqus working directory
------------------------------------------------

If you want to run Abaqus from a directory on :ref:`/fastdata <filestore>`
then you need to have the following line in your batch job submission script
just before the main ``abaqus`` command: ::

   export BAS_DISABLE_FILE_LOCKING=1

Otherwise your Abaqus job will fail and 
you will see errors like the following
in your ``my_job_name.dat`` output file: ::

    ***ERROR: An error occurred during a write access to 
              <rank=0,arg_name=outdir>my_user_job.stt file. Check the disk space 
              on your system.

This is a lie; Abaqus is failing to write the ``.stt`` file as it tries to use `file locking <https://en.wikipedia.org/wiki/File_locking>`__ 
which is not enabled on the ``/fastdata`` filesystem at present for performance reasons.
Setting the ``BAS_DISABLE_FILE_LOCKING`` environment variable to ``1`` is a Dassault Systems-approved workaround for this.



Licensed options
----------------

All available Abaqus licenses can be viewed using ``abaqus licensing r`` e.g. ::

   $ module load apps/abaqus/2017/binary
   $ abaqus licensing r

   Feature                         Version     #licenses    Expires      Vendor
   _______                         _________   _________    __________   ______
   abaqus_extended                 61.9         19          31-dec-2018  ABAQUSLM
   abaqus                          61.9         250         31-dec-2018  ABAQUSLM
   ams                             61.9         1           31-dec-2018  ABAQUSLM
   aqua                            61.9         250         31-dec-2018  ABAQUSLM
   cosim_acusolve                  61.9         1           31-dec-2018  ABAQUSLM
   cosim_direct                    61.9         1           31-dec-2018  ABAQUSLM
   cse                             61.9         1           31-dec-2018  ABAQUSLM
   design                          61.9         250         31-dec-2018  ABAQUSLM
   euler_lagrange                  61.9         1           31-dec-2018  ABAQUSLM
   gpgpu                           61.9         1           31-dec-2018  ABAQUSLM
   multiphysics                    61.9         1           31-dec-2018  ABAQUSLM
   parallel                        61.9         16384       31-dec-2018  ABAQUSLM
   sw_assoc_import                 61.9         1           31-dec-2018  ABAQUSLM
   catiav5_assoc_import            61.9         1           31-dec-2018  ABAQUSLM
   catiav5_import                  61.9         1           31-dec-2018  ABAQUSLM
   catiav6_assoc_import            61.9         1           31-dec-2018  ABAQUSLM
   tomee                           61.9         1           31-dec-2018  ABAQUSLM
   pydriver                        61.9         1           31-dec-2018  ABAQUSLM
   cae                             61.9         19          31-dec-2018  ABAQUSLM
   rtgateway                       61.9         19          31-dec-2018  ABAQUSLM
   gateway                         61.9         19          31-dec-2018  ABAQUSLM
   safe_ex_gui                     61.9         19          31-dec-2018  ABAQUSLM
   cfd                             61.9         250         31-dec-2018  ABAQUSLM
   explicit                        61.9         250         31-dec-2018  ABAQUSLM
   foundation                      61.9         250         31-dec-2018  ABAQUSLM
   simflow                         61.9         250         31-dec-2018  ABAQUSLM
   standard                        61.9         250         31-dec-2018  ABAQUSLM
   cse_token                       61.9         250         31-dec-2018  ABAQUSLM
   safe_ex_engine                  61.9         250         31-dec-2018  ABAQUSLM
   tosca_topo                      61.9         250         31-dec-2018  ABAQUSLM
   tosca_shape                     61.9         250         31-dec-2018  ABAQUSLM
   tosca_bead                      61.9         250         31-dec-2018  ABAQUSLM
   tosca_sizing                    61.9         250         31-dec-2018  ABAQUSLM
   tosca_int_abaqus                61.9         250         31-dec-2018  ABAQUSLM
   tosca_int_ansys                 61.9         250         31-dec-2018  ABAQUSLM
   tosca_int_nastran               61.9         250         31-dec-2018  ABAQUSLM
   tosca_adv_nonlinear             61.9         250         31-dec-2018  ABAQUSLM
   tosca_adv_durability            61.9         250         31-dec-2018  ABAQUSLM
   tosca_adv_morph                 61.9         250         31-dec-2018  ABAQUSLM
   tosca_smooth                    61.9         250         31-dec-2018  ABAQUSLM
   tosca_report                    61.9         250         31-dec-2018  ABAQUSLM
   tfluid_topo                     61.9         250         31-dec-2018  ABAQUSLM
   tfluid_smooth                   61.9         250         31-dec-2018  ABAQUSLM
   tfluid_parallel                 61.9         250         31-dec-2018  ABAQUSLM
   tfluid_int_ccmp                 61.9         250         31-dec-2018  ABAQUSLM
   tfluid_int_fluent               61.9         250         31-dec-2018  ABAQUSLM

Run ``abaqus licensing`` for usage info for the Abaqus licensing sub-command. Run ``abaqus licensing ru`` to see current licence usage.

Installation notes
------------------

Abaqus 2018 was installed using the
:download:`install_abaqus_2018.sh </sharc/software/install_scripts/apps/abaqus/2018/binary/install_abaqus_2018.sh>` script; the module
file is
:download:`/usr/local/modulefiles/apps/abaqus/2018/binary </sharc/software/modulefiles/apps/abaqus/2018/binary>`. 

Abaqus 2017-ga was installed using the
:download:`install_abaqus_2017-ga.sh </sharc/software/install_scripts/apps/abaqus/2017-ga/binary/install_abaqus_2017-ga.sh>` script; the module
file is
:download:`/usr/local/modulefiles/apps/abaqus/2017-ga/binary </sharc/software/modulefiles/apps/abaqus/2017-ga/binary>`. 

Abaqus 2017 was installed using the
:download:`install_abaqus_2017.sh </sharc/software/install_scripts/apps/abaqus/2017/binary/install_abaqus_2017.sh>` script; the module
file is
:download:`/usr/local/modulefiles/apps/abaqus/2017/binary </sharc/software/modulefiles/apps/abaqus/2017/binary>`. 

Abaqus 6.14-2 was installed using the
:download:`install_abaqus.sh </sharc/software/install_scripts/apps/abaqus/6.14-2/binary/install_abaqus.sh>` script; the module
file is
:download:`/usr/local/modulefiles/apps/abaqus/6.14-2/binary </sharc/software/modulefiles/apps/abaqus/6.14-2/binary>`. 

The binary installations were tested by launching ``abaqus cae`` and by using the above batch submission scripts.
Abaqus at present does not run on more than one node when using MPI due to password-less ssh being disabled across nodes on ShARC.
