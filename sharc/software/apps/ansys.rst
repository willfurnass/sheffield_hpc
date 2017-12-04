ANSYS
=====

.. sidebar:: ANSYS
   
   :Versions: 16.1, 17.2, 18.0, 18.2, 19.0
   :Dependencies: No prerequsite modules loaded. However, If using the User Defined Functions (UDF) will also need the following: For ANSYS Mechanical, Workbench, CFX and AutoDYN: Intel 14.0 or above; Compiler For Fluent: GCC 4.6.1 or above
   :URL: http://www.ansys.com 
   :Local URL: http://www.shef.ac.uk/cics/research/software/fluent


The Ansys suite of programs can be used to numerically simulate a large variety of structural and fluid dynamics problems found in many engineering, physics, medical, aeronautics and automotive industry applications.


Usage
-----

Ansys can be activated using the module files::

    module load apps/ansys/19.0/binary
    module load apps/ansys/18.2/binary
    module load apps/ansys/18.0/binary
    module load apps/ansys/17.2
    module load apps/ansys/16.1

The Ansys Workbench GUI executable is ``ansyswb``. ``ansyswb`` can be launched during an interactive session with X Window support (e.g. an interactive ``qrshx`` session).
The Ansys executable is ``mapdl`` and ``fluent`` is the executable for Fluent. Typing ``mapdl -g`` or ``fluent -g -help`` will display a list of usage options.

Batch jobs
----------

Fluent
^^^^^^

The following is an example batch submission script, ``my_job.sh``, to run ``fluent`` version 16.1 and which is submitted to the queue by typing ``qsub my_job.sh``::

    #!/bin/bash
    #$ -cwd
    #$ -l h_rt=00:30:00
    #$ -l rmem=2G
    #$ -pe mpi-rsh 8

    module load apps/ansys/19.0

    fluent 2d -i flnt.inp -g -t8 -sge -mpi=intel -rsh -sgepe mpi-rsh
	
The script requests 8 cores using the MPI parallel environment ``mpi-rsh`` with 
a runtime of 30 mins and 2 GB of real memory per core. 
The Fluent input file is ``flnt.inp``. 
**Note:** Please use the ``mpi-rsh`` parallel environment to run MPI parallel jobs for Fluent.

Legacy information
""""""""""""""""""

For older versions of Fluent (<18.0) there is a University-of-Sheffield-specific helper script, ``runfluent``, that simplified the processed of submitting Fluent jobs using a Fluent journal.
This is deprecated.  If you do want to use ``runfluent``, then you can learn about its usage by running e.g.

    module load apps/ansys/17.2
    runfluent

interactively from a worker node.

Ansys Mechanical
^^^^^^^^^^^^^^^^

Overview of execution modes
"""""""""""""""""""""""""""

Ansys Mechanical offers various ways of parallelising the running of a simulation.
You may need to experiment with the 

* type of parallelism
* number of CPU cores
* amount of memory per core

to discover the best balance between runtime and queuing time for your application.  
Note that more parallelism is not always faster: 
for example, you may find that increasing the number of threads and cores from 8 to 16 results in an increase in runtimes.

The types of parallelism can be classified as:

* multiple processes on one node (distributed memory parallelism (DMP));
* multiple processes distributed between nodes (DMP);
* a single process and multiple threads on one node (shared memory parallelism (SMP)); this usually performs worse than either of the two DMP cases;
* GPU acceleration, which can be significantly faster for single-core jobs.

It is recommended that you start with DMP parallelism on one node (e.g. < 16 cores) then 
look to see if increasing the number of cores to more than are available from a single node 
increases or reduces run-times.

Single-core batch jobs
""""""""""""""""""""""

To run a batch job using just **one CPU core** and 6 GB of RAM you need a submission script similar to the following: ::

   #!/bin/bash
   #$ -l rmem=6G
   #$ -l h_rt=03:00:00

   module load apps/ansys/19.0

   export PATH="${ANSYSROOT}/ansys/bin:${PATH}"
   mapdl -dir $PWD -b < simulation.txt > simulation.out -j simulation

Here:

* ``-dir`` tells Ansys to use a particular initial working directory (here we use ``$PWD``, the present working directory);
* ``-b`` tells Ansys to run in batch mode, which causes the input script to be listed in Ansys' output
* ``simulation.txt`` is our APDL (*ANSYS Parametric Design Language*) script
* ``simulation.out`` is the file the Ansys output will be written to
* ``simulation`` is the Ansys job name

Batch jobs with distributed memory parallelism on one node
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

To run a batch job using **multiple processes on one node** (with multiple CPU cores; here 4 cores with 3GB RAM each): ::

   #!/bin/bash
   #$ -pe smp 4
   #$ -l rmem=3G
   #$ -l h_rt=03:00:00

   module load apps/ansys/19.0

   export PATH="${ANSYSROOT}/ansys/bin:${PATH}"
   mapdl -dis -dir $PWD -np $NSLOTS -b < simulation.txt > simulation.out -j simulation

Here:

* ``-dis`` says we want to use Distributed Ansys
* ``-np`` is the number of threads Ansys is to use.
  By specifying ``$NSLOTS`` we can use the number of CPU cores specified to the scheduler (in the ``-pe smp`` line).

Batch jobs with shared memory parallelism
"""""""""""""""""""""""""""""""""""""""""

To run a batch job using **one process but multiple threads** (with multiple CPU cores **on one node**; here 4 cores with 3GB RAM each): ::

   #!/bin/bash
   #$ -pe smp 4
   #$ -l rmem=3G
   #$ -l h_rt=03:00:00

   module load apps/ansys/19.0

   export PATH="${ANSYSROOT}/ansys/bin:${PATH}"
   mapdl -dir $PWD -np $NSLOTS -b < simulation.txt > simulation.out -j simulation

Here:

* ``-np`` is the number of threads Ansys is to use.
  By specifying ``$NSLOTS`` we can use the number of CPU cores specified to the scheduler (in the ``-pe smp`` line).

Note that shared-memory parallelism often results in worse performance than distributed-memory parallelism, 
even when the allocated CPU cores are all on the same machine.

Batch jobs that use GPUs to accelerate computation
""""""""""""""""""""""""""""""""""""""""""""""""""

To run a batch job using **one CPU core and one GPU**: ::

   #!/bin/bash
   #$ -l rmem=6G
   #$ -l gpu=1
   #$ -l h_rt=03:00:00

   module load apps/ansys/19.0
   module load libs/CUDA/8.0.44/binary

   export PATH="${ANSYSROOT}/ansys/bin:${PATH}"
   mapdl -acc nvidia -na 1 -dir $PWD -b < simulation.txt > simulation.out -j simulation

Here:

* ``-acc`` is the type of accelerator (GPU) to be used.  This should always be ``nvidia`` on ShARC.
* ``-na`` signifies the number of GPUs to use.  This should typically match the number in the ``#$ -l gpu=X`` line above.

Batch jobs that distribute work between processes/machines
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

To run a batch job using **multiple processes that may be distributed between machines using MPI** (i.e. distributed memory parallelism; here we have 8 processes with 3GB RAM per process): ::

   #!/bin/bash
   #$ -pe mpi-rsh 8
   #$ -l rmem=3G
   #$ -l h_rt=03:00:00

   echo List of machines supplied to Ansys Mechanical
   machines="$(cat $PE_HOSTFILE | awk '{print $1 ":" $2}' | paste -sd:)"
   echo $machines

   module load apps/ansys/19.1

   export PATH="${ANSYSROOT}/ansys/bin:${PATH}"
   mapdl -dis -dir $PWD -usersh -mpi intelmpi -machines $machines -b < simulation.txt > simulation.txt -j simulation

None of the following should be changed:

* ``-dis`` says we want to use Distributed Ansys
* ``-mpi`` tells Mechanical that we want to use a particular version of MPI to distribute work between processes
* ``-usersh`` is necessary for MPI to work with the version of MPI we've chosen
* ``-machines`` tells Ansys's MPI how to distribute processes between nodes

.. warning:: At present only Ansys >= 19.1 is configured for efficient communication between nodes using MPI.

Legacy information
""""""""""""""""""

Previously it was recommended to use the ``runansys`` script to submit Ansys Mechanical jobs.
This command submitted an Ansys input file into the batch system and took a number of different parameters, according to your requirements.  ``runansys`` only works with Ansys Mechanical < 18.0.
To display information about how to use it, run the following from worker node:
    
    module load apps/ansys/17.2
    runansys
	
Installation notes
------------------

Ansys 19.0
^^^^^^^^^^
Installed using the :download:`install_ansys_190.sh </sharc/software/install_scripts/apps/ansys/19.0/binary/install_ansys_190.sh>` script; 
the module file is :download:`/usr/local/modulefiles/apps/ansys/19.0/binary </sharc/software/modulefiles/apps/ansys/19.0/binary>`.

The binary installations were tested by launching ``ansyswb`` and by using the above batch submission script. The ``mpi-rsh`` tight-integration parallel environment is required to run Ansys/Fluent using MPI due to password-less ssh being disabled across nodes on ShARC.

Ansys 18.2 
^^^^^^^^^^
Installed using the :download:`install_ansys_182.sh </sharc/software/install_scripts/apps/ansys/18.2/binary/install_ansys_182.sh>` script; 
the module file is :download:`/usr/local/modulefiles/apps/ansys/18.2/binary </sharc/software/modulefiles/apps/ansys/18.2/binary>`. 

Ansys 18.0 
^^^^^^^^^^

Installed using the :download:`install_ansys_180.sh </sharc/software/install_scripts/apps/ansys/18.0/binary/install_ansys_180.sh>` script; 
the module file is :download:`/usr/local/modulefiles/apps/ansys/18.0/binary </sharc/software/modulefiles/apps/ansys/18.0/binary>`. 

Ansys 17.2 
^^^^^^^^^^
Installed using the :download:`install_ansys.sh </sharc/software/install_scripts/apps/ansys/17.2/install_ansys.sh>` script; 
the module file is :download:`/usr/local/modulefiles/apps/ansys/17.2 </sharc/software/modulefiles/apps/ansys/17.2>`. 

Ansys 16.1 
^^^^^^^^^^
Installed using the :download:`install_ansys.sh </sharc/software/install_scripts/apps/ansys/16.1/install_ansys.sh>` script; 
the module file is :download:`/usr/local/modulefiles/apps/ansys/16.1 </sharc/software/modulefiles/apps/ansys/16.1>`.

Testing notes
-------------

The binary installations were tested by launching ``ansyswb`` and 
by using the above batch submission script. 
The ``mpi-rsh`` tight-integration parallel environment is required to run Ansys/Fluent using MPI due 
to password-less ssh being disabled across nodes on ShARC.

The Intel MPI bundled with Ansys can be tested using: ::

   #!/bin/bash
   #$ -pe mpi-rsh 2
   #$ -l rmem=200G
   # Or enough memory to force slots to be distributed between nodes
   machines="$(cat $PE_HOSTFILE | awk '{print $1 ":" $2}' | paste -sd:)"
   echo "Slot distribution:"
   echo $machines
   module load apps/ansys/18.2
   export PATH="${ANSYSROOT}/ansys/bin:${PATH}"
   mpitest182 -usersh -mpi intelmpi -machines $machines
