.. _embarrassing:

Embarrassingly Parallel
=======================

Overview
--------

:ref:`Slurm job arrays<array_jobs>` enable you to submit jobs that run multiple times with the same Slurm parameters.
Use the ``--array=`` Slurm argument to define array indices, e.g. ``--array=1-10,12-15``.
The ``$SLURM_ARRAY_TASK_ID`` environment variable provides each job with its
corresponding array index.

Example:

.. code-block:: slurm

  #!/bin/bash
  #SBATCH --array=1-10

  # Each job processes a different input file
  srun python my_script.py input_${SLURM_ARRAY_TASK_ID}

* Below are different templates that you can modify to suit your tasks.
* If you are uncertain about scaling up, contact the `IT Services' Research
  and Innovation team <mailto:research-it@sheffield.ac.uk>`_ at an early stage.

In scientific computing, it is often necessary to run the same program multiple times
with varying datasets or parameters.

When these runs do not depend on or communicate with each other,
they can be executed in parallel as separate Slurm jobs. This type
of parallelism is referred to as **embarrassingly parallel**.

Slurm provides a feature called **job arrays**, which allows users
to efficiently submit and manage multiple independent instances of
the same job script.

Array jobs enable you to manage large-scale workloads on the cluster.
In :ref:`parallel`, we explore alternatives.


Introduction
------------

Array jobs facilitate parallel computations. They are useful when
you need to run a job multiple times with only minor variations.
For example, you may need to execute 1000 jobs, each with a different
random seed, or apply the same operation across multiple datasets.
This can be accomplished with a single array job.

A Slurm job array consists of multiple jobs that share the same batch submission
script. The ``--array`` directive specifies how many times the script
should be executed, for instance:

.. code-block:: slurm

  #SBATCH --array=0-4

This command creates an array of five jobs (tasks) indexed from 0 to 4.
Each task is a duplicate of the submitted batch script, automatically
queued in Slurm. The ``SLURM_ARRAY_TASK_ID`` environment variable is
used to assign a unique identifier to each task, which can be leveraged
for handling input/output files.

.. image:: /images/Array-jobs.png
  :scale: 40%

.. admonition:: ``--array`` via the command line

   The ``--array`` option can also be specified as a command-line argument
   when using ``sbatch``. This is useful for managing job arrays without
   modifying the script.

.. important::

   Since array jobs create multiple identical job instances, it is
   crucial to understand their impact on the file system:

   - Does the script rely on libraries or environments stored in the working directory?
   - How much input data does each task require?
   - How much output data does each job generate?

   For example, launching an array job with hundreds of tasks that
   depend on a Python environment stored on shared storage may cause
   significant file system load due to repeated access to thousands of files.
   
   If you are unsure how your job will behave, seek guidance from the `IT Services' Research
   and Innovation team <mailto:research-it@sheffield.ac.uk>`_.

Your First Array Job
--------------------

.. include:: /referenceinfo/imports/hpc-examples-repo.rst

Let's see an array job in practice. Let's use the script names ``array_example.sh``

.. literalinclude:: /hpc-examples/examples/array/array_example.sh
   :language: slurm

Submitting the job script with ``sbatch ${HPC_EXAMPLES}/array/array_example.sh``
will return a message such as::

  Submitted batch job 5825026

This job ID belongs to the primary array job, which encompasses all
individual tasks in the array. Each task is also assigned a unique
array task ID.

As multiple jobs run simultaneously, each requires a unique output file
to prevent overwriting. By default, Slurm names the output files as
``slurm-${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.out``.

You can override this with ``--output=FILENAME``, using placeholders
``%A`` for the job ID and ``%a`` for the array task ID.

Once the jobs complete, output files will appear in your working directory:

.. code-block:: console

   $ ls
   array_example_5825026_0.out   array_example_5825026_12.out  array_example_5825026_15.out
   array_example_5825026_3.out  array_example_5825026_6.out  array_example_5825026_9.out
   array_example_5825026_10.out  array_example_5825026_13.out  array_example_5825026_1.out
   array_example_5825026_4.out  array_example_5825026_7.out  array_example.sh
   array_example_5825026_11.out  array_example_5825026_14.out  array_example_5825026_2.out
   array_example_5825026_5.out  array_example_5825026_8.out

You can inspect any output file using ``cat``::

   $ cat array_example_5825026_11.out
   I am array task number 11

.. important::

   Array indices do not need to be sequential. If specific tasks fail,
   you can re-run only those with ``--array=1,4``.
   The ``--array`` argument can also be supplied directly to ``sbatch``
   from the command line.

More Examples
-------------

The following examples demonstrate how to effectively use job arrays and leverage the ``$SLURM_ARRAY_TASK_ID`` environment variable.

- You need a clear mapping between job indices and configurations, which could be filenames, pre-defined parameter sets, or external configuration files.
- Ensure the mapping remains consistent throughout multiple job runs to avoid inconsistencies.

Processing Multiple Input Files
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Often, computations require processing different input files. The ``$SLURM_ARRAY_TASK_ID`` variable can dynamically assign files to jobs:

.. code-block:: slurm

    #!/bin/bash
    #SBATCH --time=01:00:00
    #SBATCH --mem=1G
    #SBATCH --array=0-29

    srun ./my_application -input input_data_${SLURM_ARRAY_TASK_ID}

Hardcoding Arguments in the Batch Script
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

You can explicitly set arguments within the batch script. Suppose you want to run a pi estimation simulation with five different seed values, each executing 2.5 million iterations:

Case-Based Argument Selection
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

We will use ``${HPC_EXAMPLES}/examples/array/pi_array_hardcoded_case.sh`` :

.. literalinclude:: /hpc-examples/examples/array/pi_array_hardcoded_case.sh
   :language: slurm

Submit the script with:

.. code-block:: console

   $ module load hpc-examples
   $ sbatch ${HPC_EXAMPLES}/array/pi_array_hardcoded_case.sh
   Submitted batch job 5825718

Each task produces its own output, such as:

.. code-block:: console

   $ cat pi_18.json
   {"pi_estimate": 3.1411456, "iterations": 2500000, "successes": 1963216}

Using Bash Arrays for Parameter Selection
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

An alternative approach is using Bash arrays ``${HPC_EXAMPLES}/hpc-examples/examples/array/pi_array_hardcoded_array.sh`` :

.. literalinclude:: /hpc-examples/examples/array/pi_array_hardcoded_array.sh
   :language: slurm
   

Submit the job with:

.. code-block:: console

   $ module load hpc-examples
   $ sbatch ${HPC_EXAMPLES}/array/pi_array_hardcoded_array.sh

Reading Parameters from a File
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Rather than hardcoding values, you can store them in a file and read them dynamically. For example, running ``pi.py`` with different iteration values:

Create a file named ``iterations.txt`` containing:

.. code-block:: console

   100
   1000
   50000
   1000000

We modify the script to read values using ``sed`` (see `sed <https://en.wikipedia.org/wiki/Sed>`_  and ``man sed``)
``${HPC_EXAMPLES}/array/pi_array_parameter.sh`` :

.. literalinclude:: /hpc-examples/examples/array/pi_array_parameter.sh
   :language: slurm

This approach can be extended to read multiple parameters from CSV files or similar structured data formats.

Two-Dimensional Array Scanning (Advanced) 
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

If your tasks are very short (a few minutes), launching numerous individual jobs
can lead to scheduling inefficiencies and an overwhelming number of output files.
In such cases, it is beneficial to group multiple tasks within a single array job.

.. important::

     Ideally, each array job should run for at least 30 minutes. If your tasks are
     shorter than this, consider combining multiple runs into a single job to
     reduce scheduling overhead and improve efficiency.

A simple way to achieve this is by introducing a loop inside your Slurm script.
For example, if you need to run a simulation with 50 different seed values,
you can process them in groups of 10, reducing the number of array jobs to just 5.
This significantly decreases the load on the scheduler.

An example implementation is provided in the `hpc-examples`_ script:
``examples/array/pi_array_grouped.sh``

.. literalinclude:: /hpc-examples/examples/array/pi_array_grouped.sh
  :language: slurm

Exercises
---------


Array Job Exercises
~~~~~~~~~~~~~~~~~~~

.. admonition:: Array-1: Compute n-grams with Array Jobs

      Computing n-grams across the Gutenberg-Fiction dataset can take considerable time.
      Using array jobs is an efficient way to parallelise the process. Follow along
      with this example:
      
      The following batch script ``${HPC_EXAMPLES}/ngrams/array.sh`` calculates 3-grams in 20 batches, saving each result
      to a separate file:
      
      .. literalinclude:: /hpc-examples/examples/ngrams/array.sh
         :language: slurm

      The final output now contains all computed n-grams:
      
      .. code-block:: console

            $ head -5 ngrams3-words-all.out
            30224 ["i", "don", "t"]
            18737 ["one", "of", "the"]
            15954 ["out", "of", "the"]
            14749 ["there", "was", "a"]
            13122 ["it", "was", "a"]


**Further Exercises**

.. admonition:: Array-2: Array Jobs and Random Seeds

      Create an array job that runs ``${HPC_EXAMPLES}/slurm/pi.py`` with different combinations of
      iteration counts and seed values. Save the results to separate files and keep
      the standard output (``#SBATCH --output=FILE``) distinct from the standard error (``#SBATCH --error=FILE``).

.. admonition:: Array-3: Merging Outputs

      Use the script ``${HPC_EXAMPLES}/slurm/pi_aggregation.py`` to aggregate
      results from multiple output files. This will compute a more precise estimate
      of Pi.

.. admonition:: Array-4: Applying Array Jobs to Your Own Work

      Consider your typical workload. How could you divide it into smaller, independent
      tasks that can be processed in parallel using array jobs? Would it be more
      efficient to break larger tasks into multiple smaller ones?

.. admonition:: (Advanced) Array-5: Using Advanced Indexing

      Create a job array that runs every alternate index, such as 1, 3, 5, etc. The
      `slurm sbatch manual page <https://slurm.schedmd.com/sbatch.html>`_ provides helpful
      details.

      .. dropdown:: Solution
      
           You can specify a step size for the job array using a colon and number after the range. For example: ``--array=1-X:2`` 

.. admonition:: Array-6: Varying Memory Requirements

      Construct an array job that runs ``${HPC_EXAMPLES}/slurm/memory-use.py`` with five different
      memory requirements (50M, 100M, 500M, 1000M, 5000M). Request 250M of memory
      for the array job itself. Observe whether any of the jobs fail.

      Is this an appropriate use of array jobs?

      .. dropdown:: Solution

           At a minimum, the 5G job should fail. The 500M and 1G jobs also exceed their requested memory,
           but SLURM tolerates slight overuse before terminating jobs, so they may still succeed.

           This is an incorrect use of array jobs. Arrays are intended for multiple tasks with identical
           resource requirements, as each task is allocated the same resources.

See Also
--------

* If you are uncertain about scaling up, contact the `IT Services' Research
  and Innovation team <mailto:research-it@sheffield.ac.uk>`_ at an early stage.

What's Next?
------------

The next tutorial covers :ref:`shared memory parallelism <parallel_SMP>`.



.. include:: /referenceinfo/imports/attrib_AaltoSciComp.rst

.. _hpc-examples: https://github.com/rcgsheffield/hpc-examples
