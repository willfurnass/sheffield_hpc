.. include:: warning.rst 


Flame
=====

.. sidebar:: Flame
   
   :Version: 0.17.1-openmpi-1.8.8
   :Dependancies: gcc/4.8.2 openmpi/1.8.8
   :URL: https://github.com/FLAME-HPC/xparser/archive/0.17.1.tar.gz
   :Documentation: http://flame.ac.uk/docs/


The FLAME framework is an enabling tool to create agent-based models that can be run on high performance computers (HPCs). Models are created based upon a model of computation called (extended finite) state machines. By defining agent-based models in this way the FLAME framework can automatically generate simulation programs that can run models efficiently on HPCs. The behaviour model is based upon state machines which are composed of a number of states with transition functions between those states. There is a single start state and by traversing states using the transition functions the machine executes the functions until it reaches an end state. This happens to each agent/machine as one time step or iteration is completed


Usage
-----

Flame can be activated using the module file::

    module load apps/gcc/4.8.2/flame/0.17.1-openmpi-1.8.8


Note: The module file also loads libtool, libmboard, and the compilers gcc/4.8.3 & openmpi/1.8.8. Libmboard is the communication library used by Flame simulation programs.

Running Flame in parallel on multiple cores requires the batch script to be placed within the model folder. The batch script must contain the following instructions

.. code-block:: bash

   #!/bin/bash
   #$ -pe openmpi-ib no_cores
   #$ -cwd
   #$ -j y
   #$ -m bea
   #$ -M email_address
   module load apps/gcc/4.8.2/flame/0.17.1-openmpi-1.8.8
   make clean
   cp /usr/local/packages6/apps/gcc/4.8.2/flame/0.17.1-openmpi-1.8.8/xparser/0.17.1/*.tmpl .
   xparser model.xml -p
   make
   mpirun -np no_cores ./main no_its its/0.xml

where

no_cores = the number of computer cores required
no_its = the number of time steps for the simulation
0.xml = the initial state of the system/agents. This file is placed within the its folder within the model folder


Test
----

#. Download and unpack the following example models using:

   .. code-block:: bash
   
      cd some_directory
      wget http://flame.ac.uk/docs/zip/tutorial_models.zip
      unzip tutorial_models.zip
      cd tutorial_models/model_0x

#. create the ``its`` folder within the model_0* folder, and move the 0.xml file to the ``its`` folder

#. Create a batch job submission script called ``test.sge`` containing the 12 lines of code above. Replace model.xml (line 10) with the actual xml model name. The script should be in the model folder.

#. Submit the job using ``qsub test.sge``


Installation notes
------------------

Flame was compiled using the
:download:`install_flame.sh </iceberg/software/install_scripts/apps/gcc/4.8.2/flame/0.17.1-openmpi-1.8.8/install_flame.sh>` script, the module
file is
:download:`0.17.1-openmpi-1.8.8 </iceberg/software/modulefiles/apps/gcc/4.8.2/flame/0.17.1-openmpi-1.8.8>`.
