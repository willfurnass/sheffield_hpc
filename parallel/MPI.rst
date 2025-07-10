.. _parallel_MPI:

MPI for Multi-node and Parallel Jobs
====================================

Overview
--------

* Confirm that your application is built to support MPI.

* Compile your code using MPI-aware compilers. Be sure to load the same modules in your job script.

* To assign multiple tasks on a single machine, use ``--nodes=1`` and ``--ntasks=n``.

* Launch the program using ``srun`` if you're relying on system-installed MPI modules, or use ``mpirun`` for custom-built MPI.

* To distribute tasks evenly across multiple machines, combine ``--nodes=N`` with ``--ntasks-per-node=n`` to yield a total of **N x n** tasks.

* Always check resource usage with tools like ``seff JOBID`` to ensure efficient use of allocations.

* If you are uncertain about scaling up, contact the `IT Services' Research
  and Innovation team <mailto:research-it@sheffield.ac.uk>`_ at an early stage.


MPI allows a program to run concurrently across many cluster nodes, although it typically requires more specialised programming.


What is MPI?
------------

The Message Passing Interface is a standard for passing data and other messages between running `processes <https://en.wikipedia.org/wiki/Process_(computing)>`_
which may or may not be on a single computer.
It is commonly used on computer clusters as a means by which a set of related processes can work together in parallel on one or more tasks.
These strands (processes) must therefore communicate data and other information by passing messages between each other.

MPI is used on systems ranging from a few interconnected `Raspberry Pi's <http://thenewstack.io/installing-mpi-python-raspberry-pi-cluster-runs-docker/>`_ through to
the UK's national supercomputer, `Archer <http://www.archer2.ac.uk/>`_.

.. _mpi_impl:

MPI Implementations
-------------------
The `Message Passing Interface (MPI) <http://mpi-forum.org/>`_ itself is just a *specification* for a message passing library.

There are multiple implementations of this specification, each produced by a different organisation,
including `OpenMPI <https://www.open-mpi.org/>`_ and `Intel MPI <https://www.intel.com/content/www/us/en/developer/tools/oneapi/mpi-library.html>`_.
This documentation includes information on the MPI implementations available on :ref:`Stanage <stanage-parallel>` and :ref:`Bessemer <bessemer-parallel>`.
On the Stanage cluster these implementations have been compiled in a way that allows them to make optimal use of the high-speed network infrastructure (OmniPath).
If you are not sure which implementation to use then try the latest available version of OpenMPI.

Batch MPI
---------
To use MPI you need use ``module load`` to activate a particular :ref:`MPI implementation <mpi_impl>`
(or ``module load`` an application that itself loads an MPI implementation behind the scenes).

Here is an example that requests 4 *slots* (CPU cores) with 8GB of RAM per slot then runs a program called ``executable``
in the current directory using the OpenMPI library (version 4.1.4, built using version 12.2.0 of the gcc compiler).
It is assumed that ``executable`` was previously compiled using that exact same MPI library.

.. code-block:: console

   #!/bin/bash
   # Request one node
   #SBATCH --nodes=1
   # Request 4 cores per node
   #SBATCH ntasks=4
   # Request 8GB of RAM per node
   #SBATCH mem=8G

   # Load a MPI library
   module load OpenMPI/4.1.4-GCC-12.2.0

   # Run a program previously compiled using that specific MPI library
   srun --export=ALL ./executable

Unlike shared memory models, MPI requires programs to explicitly send and receive data between tasks.
Most applications must be written with MPI support from the outset,
so standard serial code won’t benefit from MPI unless rewritten accordingly.

.. image:: /images/MPI.png
  :scale: 40%

MPI programs usually follow this pattern:

1. The same executable is launched in several separate processes.
2. All processes connect to one another through MPI.
3. Each process is assigned a unique identifier called a "rank".
4. Each rank performs a specific portion of the work. Rank 0 often handles I/O and status messages.
5. The MPI environment is closed after execution ends.

If you're using MPI modules provided by the system, Slurm communicates rank and task details using a
library called `PMIx <https://pmix.org/>`_. This may not work seamlessly with external MPI builds.

Building and Executing MPI Programs
-----------------------------------

Compiling with MPI
~~~~~~~~~~~~~~~~~~

Choose an MPI implementation for building your application. Several are available, all conforming to the MPI standard.
We suggest using the latest version of OpenMPI for compatibility with the cluster environment.
For information on other installed versions, see :ref:`stanage-parallel`.


Requesting MPI Resources in Slurm
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

To allocate resources for an MPI job, the typical format is: ``--nodes=1 --ntasks=N``.
This ensures all MPI tasks run on a single machine—ideal for communication-intensive programs.

When scaling to multiple nodes, use: ``--nodes=N --ntasks-per-node=n``. This launches **N x n** tasks,
balancing them across machines. Each task gets 1 CPU by default. To increase this (if your programme supports it), see the section on
:ref:`hybrid parallel models <mpi-hybrid-parallelisation>`.

.. _pi-mpi-example:

Example: MPI Program to Estimate Pi
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

We'll use the ``pi-mpi.c`` example, which estimates π using a Monte Carlo method, and supports multiple MPI tasks.

Start by compiling with the MPI module:

.. code-block:: console

    module load hpc-examples
    module load OpenMPI
    mpicc -o pi-mpi ${HPC_EXAMPLES}/slurm/pi-mpi.c

Run interactively:

.. code-block:: bash

    srun --nodes=1 --ntasks=2 --time=00:10:00 --mem=500M ./pi-mpi 1000000

Or via a Slurm script (``pi-mpi.sh``):

.. code-block:: slurm

    #!/bin/bash
    #SBATCH --time=00:10:00
    #SBATCH --mem=500M
    #SBATCH --output=pi-mpi.out
    #SBATCH --nodes=1
    #SBATCH --ntasks=2

    module load OpenMPI
    srun --export=ALL ./pi-mpi 1000000

Submit with:

.. code-block:: slurm

    sbatch pi-mpi.sh

You can inspect the output file using ``cat``:

.. code-block:: console

  $ cat pi-mpi.out
  node032.pri.stanage.alces.network: This is rank 1 doing 500000 trials
  Calculating pi using 1000000 stochastic trials
  node032.pri.stanage.alces.network: This is rank 0 doing 500000 trials
  Throws: 785491 / 1000000 Pi: 3.141964

.. important::

   Here we didn’t specify an OpenMPI version, so the system default was used.
   However, for reproducibility and to avoid runtime errors, **always load the same version of OpenMPI**
   as you used when compiling the program.
   Mismatched major versions (e.g. 3.x vs 4.x) can cause MPI initialisation errors or crashes,
   and even minor differences can affect runtime behaviour.

Special Cases and Tips
-----------------------

Ranks Not Detected
~~~~~~~~~~~~~~~~~~

When using MPI libraries outside the system default, ranks may not be recognised automatically. Add the following to your job script if needed:

.. code-block:: slurm

    export SLURM_MPI_TYPE=pmix_v2

*Note: We suggest* ``pmix_v2`` *here for broad compatibility with different MPI builds, but newer MPI libraries may also support* ``pmix_v4``.

Performance Monitoring
~~~~~~~~~~~~~~~~~~~~~~

.. include:: /referenceinfo/imports/scheduler/SLURM/common_commands/seff_usage_import.rst

.. _mpi-hybrid-parallelisation:

Hybrid Parallelism: MPI + Shared Memory
---------------------------------------

Some applications benefit from both task-based and thread-based parallelism (e.g., MPI + OpenMP). This is often called hybrid parallelism.

.. warning::

   Hybrid parallelism is generally only beneficial on very large-scale resources (thousands to tens of thousands of cores).
   For smaller jobs, the overhead often outweighs any performance gains.

To request multiple threads per MPI task, use: ``--ntasks-per-node=n --cpus-per-task=C``

Ensure the total cores requested per node (**n × C**) does not exceed the node capacity.

Exercises
---------

.. admonition:: Exercise 1: Understanding Basic MPI Options

    Try these commands and observe their behaviour:

    .. code-block:: console

        srun --cpus-per-task=4 hostname
        srun --ntasks=4 hostname
        srun --nodes=4 hostname

    .. dropdown:: Solution

       - The first uses 4 CPUs in one task.
       - The second starts 4 tasks, one CPU each.
       - The third won't work, as Stanage doesn't permit interactive jobs spanning multiple nodes.
         Therefore, you'll need to submit a batch script ``submit.sh`` to Slurm:

       .. code-block:: slurm

          #!/bin/bash
          #SBATCH --nodes=4
          #SBATCH --ntasks-per-node=1
          #SBATCH --time=00:01:00
          #SBATCH --mem=2G
          #SBATCH --output=4node-test.out

          srun hostname

       To submit: ``sbatch submit.sh``

       The output file will show that the tasks were spread across 4 nodes.

.. admonition:: Exercise 2: Run MPI with Various Slurm Options

    Try the ``pi-mpi.c`` :ref:`example <pi-mpi-example>` using:

    1. ``--ntasks=4``
    2. ``--ntasks-per-node=4``
    3. ``--nodes=2 --ntasks-per-node=2``

    .. dropdown:: Solution

      You can test with:

      .. code-block:: bash

         module load OpenMPI
         srun --export=ALL --ntasks=4 --time=00:10:00 --mem=500M ./pi-mpi 2000000000
         srun --export=ALL --ntasks-per-node=4 --time=00:10:00 --mem=500M ./pi-mpi 2000000000

      Again, we will need to submit a script to Slurm for a multi-node job:

      .. code-block:: slurm

         #!/bin/bash
         #SBATCH --nodes=2
         #SBATCH --ntasks-per-node=2
         #SBATCH --time=00:10:00
         #SBATCH --mem=500M
         #SBATCH --output=pi-mpi-multi-node-test.out

         module load OpenMPI
         srun --export=ALL ./pi-mpi 2000000000

      The difference is how tasks are distributed across nodes.

      - In the first example, all 4 tasks run on a single node (default behaviour).
      - In the second example, Slurm packs 4 tasks onto a single node explicitly.
      - In the third example, Slurm distributes 2 tasks to each of two nodes.

      These examples help you see how **task placement** influences performance.

      You can inspect the efficiency of these jobs with ``seff JOBID``:

      .. csv-table::
         :header-rows: 1
         :delim: |

         Distribution                   | Wall Time  | CPU Time  | CPU Efficiency (%)
         ``ntasks=4``                   | 00:08      | 00:32     | 96.88
         ``ntasks-per-node=4``          | 00:09      | 00:36     | 86.11
         ``ntasks=2 ntasks-per-node=2`` | 00:17      | 01:08     | 60.29

      - The slight difference between the first two cases is likely just due to small startup and timing artefacts, typical for short jobs.
      - The significant drop in efficiency in the third case is expected — it reflects the overhead of communicating across nodes, especially for small jobs.


      .. tip::

         - Multi-node jobs have **more communication overhead** due to network latency.
         - For larger jobs, consider scaling tests to find the “sweet spot” between cores, memory, and nodes.

.. admonition:: Exercise 3: Can Your Code Use MPI?

    Look at your code’s documentation or output. Keywords that hint at MPI support include:

    * MPI
    * mpirun
    * mpiexec
    * distributed
    * rank

.. admonition:: Exercise 4:

  Explore our documentation pages on parallel implementations on  :ref:`Stanage <stanage-parallel>` and :ref:`Bessemer <bessemer-parallel>`.
  Pay attention to examples, best practices, and cluster-specific tweaks -
  they’ll give you a head start in deploying MPI effectively on our systems.

MPI Training
------------
Training courses from the national supercomputing centre are available `here <https://www.archer2.ac.uk/training/courses/210000-mpi-self-service/>`_

What’s Next?
------------

The next guide introduces :ref:`GPU parallelism <GPU_computing>` and how to use GPUs on the cluster.


.. include:: /referenceinfo/imports/attrib_AaltoSciComp.rst
