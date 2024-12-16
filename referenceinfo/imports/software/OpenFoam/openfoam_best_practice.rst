

By default, each process in an OpenFOAM parallel simulation writes one file for each output field at each output time.

.. code-block:: console

    number of files = number of output fields x number of output times x number of processes

This can lead to a large number of small files. For example, a job with 256 cores solving 20 variables and saving 100 time steps will end up with 512,000 result files.

As OpenFOAM interacts with each one of these files, it triggers I/O operations such as reads or writes. Having this many I/O operations happening in parallel can lead to a number of performance and stability problems. Below are some of the main issues that might arise when dealing with so many I/O operations:

#. **Resource Contention**    
    * Disk/Network Bandwidth Saturation: Parallel I/O operations may overwhelm the available bandwidth of storage devices or networks. This leads to congestion and increases the time to complete each I/O operation.

#. **Increased Latency** 
    * I/O Queue Overload: When too many I/O requests are submitted simultaneously, they often get queued. This can increase the latency of each operation as they have to wait longer to be serviced.

#. **System Instability** 
    * Memory Pressure: Each I/O operation consumes memory (buffers, caches, control blocks). Too many parallel operations can exhaust system memory, leading to excessive swapping or out-of-memory (OOM) conditions, causing crashes. In HPC Slurm jobs, this may result in out-of-memory failures.

#. **Starvation of Critical Tasks**
    * Important Tasks Delayed: With too many simultaneous operations, more critical I/O tasks (like responding to user requests or handling time-sensitive operations) might be delayed or starved of resources, degrading system responsiveness.

#. **Deadlocks and Race Conditions**
    * Complexity in Synchronization: The more I/O operations you run in parallel, the more complex the synchronization mechanisms (locks, semaphores, etc.) become. This increases the risk of deadlocks or race conditions, leading to bugs and system instability.

#. **Increased Error Rates**
    * Timeouts and Retries: With too many operations in progress, some may time out due to resource starvation or congestion, leading to increased retries and error rates. This can further add to the system load, worsening the situation.

#. **Diminishing Returns**
    * Limited Parallelism Gains: Depending on the system architecture (e.g., SSDs vs. HDDs, network vs. local I/O), the benefit of adding more parallel I/O operations decreases after a certain point. This results in diminishing returns in performance improvements, and in some cases, can degrade performance.

#. **Poor User Experience**
    * Slow Application Performance: Excessive parallel I/O operations can slow down the entire system, resulting in sluggish application performance, unresponsiveness, and increased wait times for users.

#. **Heat and Power Consumption**
    * Increased System Load: Constant high I/O activity can cause components like CPUs, memory, and storage devices to run at high loads for extended periods, which may increase power consumption and generate more heat, potentially leading to hardware wear or thermal throttling.

Some of these issues can be remedied by requesting more resources, such as additional run time and memory, but this approach is wasteful, inefficient, and will likely result in longer job queues.

However, some of these issues may have more severe consequences, such as slowing down other users' jobs in our multi-tenant HPC system or even node failures.

These issues are magnified on :ref:`Fastdata areas (Lustre) <fastdata_dir>` which are designed for large files and handle small files poorly.

**Solutions And Best Practice**

As detailed above, OpenFOAM's default settings create many issues for the HPC system. Below are some best practices and solutions to address these problems:

**Collated file handler**

The concept of "collating" OpenFOAM output is to merge the results of a variable corresponding to a group of processors (ranks) into a single file. For example, if the "collating" groups are formed every 32 ranks, then results from processors 0 to 31 are merged together, and results from processors 32 to 63 are merged together, and so on. (This feature is only available since OpenFOAM-org-6 and OpenFOAM-v1806 versions.)

.. code-block:: console

    number of files = (number of output fields x number of output times x number of processes) / group_count_variable

    
**Settings part-1: Define the file handler with -fileHandler collated**

The collated file handler can be indicated in the ``<case>/system/controlDict`` file:

.. code-block:: console

    //...
    OptimisationSwitches
    {
        fileHandler collated;
    }
    //...

Or via the command line when executing all solvers/tools:

.. code-block:: console

    $ decomposePar -fileHandler collated

**Settings part-2: Define the size of the groups of ranks to be collated with FOAM_IORANKS**

Use the following syntax to define groups of ranks of size G:

.. code-block:: console

    $ export FOAM_IORANKS='(0 G 2G ... mG)'

For example, for grouping every 16 ranks in a case with 64 processors, use:

.. code-block:: console

    $ export FOAM_IORANKS='(0 16 32 48 64)'

Example of the resulting directory structure:

.. code-block:: console

    $ export FOAM_IORANKS='(0 16 32 48 64)'

    $ decomposePar  -fileHandler collated

    ....

    $ ls
    0       constant                processors64_32-63      system
    0.orig  processors64_0-15        processors64_16-31      Allrun  

**Best practices**

If you find that reading and writing files takes up a significant fraction of your jobs time, you can modify the input and/or output settings in ``controlDict``. Some of these suggestions may not be feasible depending on your analysis needs:

1.  **Increase writeInterval**:
    This is a ``controlDict`` parameter that controls how often time directories are stored. A lower value means fewer result files are kept on disk. For steady-state problems, only one directory is necessary, and previous time step data can be overwritten (purgeWrite 1 is recommended).

2.  **Use binary format for fields**: ``writeFormat binary``:
    Writing output in binary format is faster than ASCII. Although ASCII is human-readable, it affects performance with large simulations. Set this option in ``controlDict``.

3.  **For steady-state solutions, overwrite output at each time**: ``purgeWrite 1``.

4.  **Disable reading dictionaries at every time step**: ``runTimeModifiable no``.
    This prevents unnecessary overhead, improving performance by turning off runtime modification of parameters.
