.. _sched_batch:

Batch jobs
==========

The power of the clusters really comes from the *batch job* queue submission process.
Using this system, you write a script which
requests various resources,
initializes the computational environment and then
executes your program(s).
The scheduler will run your job when resources are available.
As the task is running, the terminal output and any errors are captured and
saved to disk, so that you can see the output and verify the execution of the
task.

Any task that can be executed without any user intervention while it is running
can be submitted as a batch job.
This typically excludes jobs that require a Graphical User Interface (GUI).
However, many common GUI applications such as Ansys or MATLAB can also be used without their GUIs.

When you submit a batch job,
you provide an executable file that will be run by the scheduler.
This is normally a bash script file which provides commands and options to the program you are using.
Once you have a script file, or other executable file, you can submit it to the queue by running:

* **SGE** ::

    qsub myscript.sh

* **Slurm** ::

    sbatch myscript.sh

You should then be informed of the unique identifier of your job e.g.: ::

    Submitted batch job 1226

When your job gets to the front of a job queue and starts running textual output from the job will, by default, be written to one or more files in the directory you submitted the job from.

So what goes in ``myscript.sh``?
Here is an example SGE batch submission script that runs a fictitious program called ``foo``:

   .. code-block:: bash

    #!/bin/bash
    # Request 5 gigabytes of real memory (mem)
    #$ -l rmem=5G

    # load the module for the program we want to run
    module load apps/gcc/foo

    # Run the program foo with input foo.dat
    # and output foo.res
    foo < foo.dat > foo.res

To use Slurm the equivalent batch submission script would be:

   .. code-block:: bash

    #!/bin/bash
    # Request 5 gigabytes of real memory (mem)
    #SBATCH --mem=5G

    # load the module for the program we want to run
    module load apps/gcc/foo

    # Run the program foo with input foo.dat
    # and output foo.res
    foo < foo.dat > foo.res


Some things to note:

* The first line always needs to be ``#!`` then the full path to a program that the scheduler will run to interpret the rest of the script. Most job scripts are (bash) shell scripts so this line is typically ``#!/bin/bash``
* Comments start with a ``#`` (if the first line is ``#!/bin/bash``)
* :ref:`Scheduler options <sched_options>`, such as the amount of memory requested, start with:
  * ``#$`` when using SGE
  * ``#SBATCH`` when using Slurm
* You will often require one or more ``module`` commands in your submission file. 
  These :ref:`make programs and libraries available to your scripts <env_modules>`
  Many applications and libraries are available as modules on 
  :ref:`ShARC <sharc-software>`, :ref:`Bessemer <bessemer-software>` and :ref:`iceberg <iceberg-software>`.

Here is a more complex example that requests more resources:

Using **SGE:**

   .. code-block:: bash

    #!/bin/bash
    # Request 16 gigabytes of real memory (RAM)
    #$ -l rmem=16G
    # Request 4 cores in an OpenMP environment
    #$ -pe openmp 4
    # Email notifications to me@somedomain.com
    #$ -M me@somedomain.com
    # Email notifications if the job aborts
    #$ -m a

    # Load the modules required by our program
    module load compilers/gcc/5.2
    module load apps/gcc/foo

    # Set the OPENMP_NUM_THREADS environment variable to 4
    export OMP_NUM_THREADS=4

    # Run the program foo with input foo.dat
    # and output foo.res
    foo < foo.dat > foo.res

Using **Slurm:**

   .. code-block:: bash

    #!/bin/bash
    # Request 16 gigabytes of real memory (RAM)
    #SBATCH --mem=16G
    # Request 4 cores 
    #SBATCH -c 4
    # Email notifications to me@somedomain.com
    #SBATCH --mail-user=me@somedomain.com
    # Email notifications if the job fails
    #SBATCH --mail-type=FAIL

    # Load the modules required by our program
    module load compilers/gcc/5.2
    module load apps/gcc/foo

    # Set the OPENMP_NUM_THREADS environment variable to 4
    export OMP_NUM_THREADS=4

    # Run the program foo with input foo.dat
    # and output foo.res
    foo < foo.dat > foo.res
