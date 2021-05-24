.. _ansys-sharc-fluent:

.. include:: ../ansys/sharc-sidebar.rst

Fluent
========================

.. contents::
    :depth: 3

----------------

ANSYS Fluent is a Computational Fluid Dynamics (CFD) code for modelling fluid flow, heat transfer, mass transfer and chemical reactions.
Fluent has interfaces to other pre and post processing packages supplied by ANSYS such as ICEMCFD or Ensight and can incorporate user developed models via user defined functions.
Fluent can make use of  built in :ref:`MPI <parallel_MPI>` to utilize multiple cross node CPU and can scale to hundreds of cores.

----------------

.. include:: ../ansys/module-load-list.rst

--------------------

Interactive jobs
----------------

While using a X11 GUI forwarding supported SSH client, an interactive session can be started on ShARC with the ``qrshx`` command which supports graphical applications.
You can load an ANSYS module above and then start the fluent program by running the ``fluent`` command.

If desired, the ANSYS Workbench GUI executable can be launched with the  ``ansyswb`` command.
To use more than a single core, you should write a batch job script and fluent journal file for submission to the batch queues.


.. caution::

  Warning for ANSYS versions >= 18.0 using the command ``fluent`` results in an unresponsive fluent launcher. To launch fluent and bypass the launcher use ``fluent dim`` where dim = 2d, 2ddp, 3d or 3ddp or unset the following environment variables before running the command::

    unset SGE_TASK_ID
    unset RESTARTED

--------------------

Batch jobs
----------

To submit Fluent jobs with larger resources a batch job must be submitted to the scheduler. A Fluent batch job is composed of a Fluent journal file which tells Fluent what to do and submission script which tells the scheduler how much resource to request and how to start Fluent.

Fluent Journal files
^^^^^^^^^^^^^^^^^^^^

A Fluent journal file is effectively a sequence of TUI and/or Scheme commands that you’ll supply to Fluent instead of using the GUI / TUI.
A more thorough explanation anmd tutorial on how to make a Fluent journal file can be found on the following page:  :ref:`Writing Fluent journal files. <writing-fluent-journal-files>`

.. important::

  It is critical that you consult with the :ref:`Fluent journal syntax <writing-fluent-journal-files>` documentation linked above in order to ensure features such as autosaving and intermediate save states are turned on and functioning correctly!

Batch Submission Script
^^^^^^^^^^^^^^^^^^^^^^^

Fluent is capable of running in both :ref:`MPI <parallel_MPI>` and :ref:`SMP <parallel_SMP>` parallel environments but will use its in-build MPI communications for both.
On ShARC, cross process communication must use the RSH protocol instead of SSH.
This necessitates the use of either the ``smp`` (up to 16 cores on a single node only) or ``mpi-rsh`` (as many cores as desired across many nodes) parallel processing environments.


Sample MPI Fluent Scheduler Job Script
"""""""""""""""""""""""""""""""""""""""""

The following is an example batch submission script, ``cfd_job.sh``, to run the executable ``fluent`` with input journal file ``test.jou``, and carry out a 2D double precision CFD simulation.
The script requests 8 cores using the MPI parallel environment ``mpi-rsh`` with a runtime of 30 mins and 2 GB of real memory per core. The Fluent input journal file is ``test.jou``.

.. hint::

    * Use of the ``#$ -V`` SGE option will instruct SGE to import your current terminal environment variables to be imported - **CAUTION** - this may not be desirable.
    * Use of the ``mpi-rsh`` parallel environment to run MPI parallel jobs for Ansys is required if using more than 16 cores on ShARC.
    * **$NSLOTS** is a Sun of Grid Engine variable which will return the requested number of cores.
    * **-mpi=intel** instructs Fluent to use the in-built Intel MPI communications - important for using the high performance :ref:`Omnipath networking <sharc-network-specs>` between nodes.
    * **-rsh** tells Fluent to use RSH instead of SSH.
    * **-sge** forces Fluent to recognise job submission via SGE.
    * **-sgepe** selects the *mpi-rsh* SGE parallel environment.
    * **-pib.infinipath** instructs Fluent to use the high performance Omnipath networking. :ref:`Omnipath networking <sharc-network-specs>`
    * **-gu** and **-driver null** instructs Fluent that it will be running with no GUI to avoid errors caused by plot / figure export.

.. code-block:: bash

    #!/bin/bash
    #$ -V
    #$ -cwd
    #$ -M joe.bloggs@sheffield.ac.uk
    #$ -m abe
    #$ -l h_rt=00:30:00
    #$ -l rmem=2G
    #$ -pe mpi-rsh 8
    #$ -N JobName

    module load apps/ansys/20.2/binary

    fluent 2ddp -i test.jou -gu -t$NSLOTS -mpi=intel -rsh -sgepe mpi-rsh -sge -pib.infinipath -driver null

-----------------------

Sample SMP Fluent Scheduler Job Script
""""""""""""""""""""""""""""""""""""""""

The following is an example batch submission script, ``cfd_job.sh``, to run the executable ``fluent`` with input journal file ``test.jou``, and carry out a 2D double precision CFD simulation. The script requests 8 cores using the SMP parallel environment ``smp`` with a runtime of 30 mins and 2 GB of real memory per core. The Fluent input journal file is ``test.jou``. The infinipath option is omitted as this job runs within a single node.

.. code-block:: bash

    #!/bin/bash
    #$ -V
    #$ -cwd
    #$ -M joe.bloggs@sheffield.ac.uk
    #$ -m abe
    #$ -l h_rt=00:30:00
    #$ -l rmem=2G
    #$ -pe smp 8
    #$ -N JobName

    module load apps/ansys/20.2/binary

    fluent 2ddp -i test.jou -g -t$NSLOTS -mpi=intel -rsh -sgepe smp -sge -driver null


Further details about how to construct batch jobs can be found on the :ref:`batch submission guide <submit-batch>` page

The job is submitted to the queue by typing:

.. code-block:: bash

    qsub cfd_job.sh

----------------

.. include:: ../../../../referenceinfo/ANSYS/fluent/export-fluent-plots-while-using-batch-jobs.rst
