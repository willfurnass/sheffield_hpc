.. include:: warning.rst 

.. _abaqus_iceberg:

Abaqus
======

.. sidebar:: Abaqus

   :Versions:  2017, 6.14, 6.13, 6.12 and 6.11
   :Support Level: FULL
   :Dependancies: Intel Compiler (for user subroutines)
   :URL: http://www.3ds.com/products-services/simulia/products/abaqus/
   :Local URL: https://www.sheffield.ac.uk/it-services/research/software/abaqus

Abaqus is a software suite for Finite Element Analysis (FEA) developed by Dassault Systèmes.


Interactive usage
-----------------

After connecting to iceberg (see :ref:`ssh`),  start an interactive session with the :code:`qsh` command. 
Alternatively, if you require more memory, for example 16 gigabytes, use the command :code:`qsh -l rmem=16G` 

Make a specific version available with one of the following commands:

.. code-block:: none

      module load apps/abaqus/2017/binary
      module load apps/abaqus/614
      module load apps/abaqus/613
      module load apps/abaqus/612
      module load apps/abaqus/611

Command-line access
^^^^^^^^^^^^^^^^^^^

To access Abaqus' command-line interface:

.. code-block:: sh

    abaqus

Graphical interface
^^^^^^^^^^^^^^^^^^^

The recommended way to start Abaqus' **graphical interface** depends on 

* which version of Abaqus has been loaded
* whether you are using a :ref:`(hardware-accelerated) 'qsh-vis' session <iceberg-hw-accel-gfx>`
  which can **greatly improve graphical performance** but requires a graphics card, which are typically in demand.

Without qsh-vis (> 6.13)
""""""""""""""""""""""""

.. code-block:: sh

    abaqus cae -mesa

With qsh-vis (<= 2017)
""""""""""""""""""""""

.. code-block:: sh

    abaqus cae

With qsh-vis (> 6.13 and < 2017)
""""""""""""""""""""""""""""""""

.. code-block:: sh

    $ABAQUSCOMMAND cae

With qsh-vis (6.13)
"""""""""""""""""""

.. code-block:: sh

    abq6133 cae


Example problems
----------------
Abaqus contains a large number of example problems which can be used to become familiar with Abaqus on the system. 
These example problems are described in the Abaqus documentation, 
and can be obtained using the Abaqus ``fetch`` command. 
For example, after loading the Abaqus module enter the following at the command line to 
extract the input file for test problem ``s4d``:

.. code-block:: sh

    abaqus fetch job=s4d

This will extract the input file ``s4d.inp``. 
To run the computation defined by this input file replace ``input=myabaqusjob`` with ``input=s4d`` in the commands and scripts below.

Batch jobs 
----------

Single-core job
^^^^^^^^^^^^^^^

In this example, we will run the ``s4d.inp`` file on a single core using 8 Gigabytes of memory.  
After connecting to iceberg (see :ref:`ssh`), 
start an interactive sesssion with the :code:`qrsh` command.

Load version 2017 of Abaqus and fetch the ``s4d`` example by running the following commands:

.. code-block:: sh

    module load apps/abaqus/2017/binary
    abaqus fetch job=s4d

Now, you need to write a batch submission file. We assume you'll call this :code:`my_job.sge`:

.. code-block:: sh

    #!/bin/bash
    #$ -cwd
    #$ -l rmem=8G

    module load apps/abaqus/2017/binary

    abaqus job=my_job input=s4d.inp scratch=$TMPDIR memory="8gb" interactive

Submit the job with:

.. code-block:: sh

    qsub my_job.sge

* We have requested 8 gigabytes of memory in the above job. The ``memory="8gb"`` switch tells abaqus to use 8 gigabytes. 
  The ``#$ -l rmem=8G`` tells the system to reserve 8 gigabytes of real memory.  Make sure that the ``memory=`` and ``rmem=`` values match.
* Note the word ``interactive`` at the end of the abaqus command. Your job will not run without it.

Single-core job with user subroutine
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

In this example, we will fetch a simulation from Abaqus' built in set of problems that 
makes use of user subroutines (UMATs) and run it in batch on a single core.  
After connecting to iceberg (see :ref:`ssh`),  
start an interactive session with the :code:`qrsh` command.

Load version 2017 of Abaqus and fetch the ``umatmst3`` example by running the following commands:

.. code-block:: sh

    module load apps/abaqus/2017/binary
    abaqus fetch job=umatmst3*

This will produce two files: 

* The input file ``umatmst3.inp`` 
* the Fortran user subroutine ``umatmst3.f``

Now, you need to write a batch submission file. We assume you'll call this :code:`my_user_job.sge`:

.. code-block:: sh

    #!/bin/bash
    #$ -cwd
    #$ -l rmem=8G

    module load apps/abaqus/2017/binary
    module load $ABAQCOMPVER

    abaqus job=my_user_job input=umatmst3.inp user=umatmst3.f scratch=$TMPDIR memory="8gb" interactive

Submit the job with: 

.. code-block:: sh

    qsub my_user_job.sge

Important notes:

* In order to use user subroutines, it is necessary to load the module for a particular version of the :ref:`Intel compiler <iceberg_intel_compilers>`.
  The name of the module file for the most appropriate Intel compiler is stored in the ``ABAQCOMPVER`` environment variable.
* The user subroutine itself is passed to Abaqus with the switch ``user=umatmst3.f``.
* The notes for the previous single-core batch job example still apply.

Multi-core job (on a single node)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

To distribute the workload associated with an Abaqus job between say 4 CPU cores on the same worker node
we need a batch job submission script ``my_4_core_job.sge`` like the following:

.. code-block:: sh

    #!/bin/bash
    #$ -cwd
    #$ -l rmem=8G
    #$ -pe openmp 4

    module load apps/abaqus/2017/binary

    abaqus job=my_job input=s4d.inp mp_mode=threads cpus=$NSLOTS scratch=$TMPDIR memory="32gb" interactive

Again, submit the job with: 

.. code-block:: sh

    qsub my_4_core_job.sge

Important notes:

* We specify the **number of CPU cores** using ``-pe openmp 4`` near the top of the script.
  We tell Abaqus to distribute the work using ``cpus=$NSLOTS`` where 
  ``NSLOTS`` is a variable automatically set by the job scheduler to be 
  the same as the number at the end of the ``-pe openmp`` line.
* Here we request a job with 8GB of real **memory per CPU core** (``-l rmem=8G``)
  but Abaqus itself needs to be told the **total amount of memory available** (``memory="32gb"``)
* The notes for the previous single-core batch job example still apply.

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