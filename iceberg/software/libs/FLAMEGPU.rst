.. Warning:: 
    Iceberg reaches end-of-life on **30th November 2020.**
    If you are running jobs on Iceberg then you need to take urgent action to ensure that your jobs/scripts will run on ShARC or Bessemer. 
 
    If you have never used ShARC or Bessemer then now is the time to test your scripts.
    Not all software on Iceberg is available on ShARC/Bessemer. 

.. _flamegpu_iceberg:

.. highlight:: bash

FLAMEGPU
========

.. sidebar:: FLAME GPU

   :URL: http://www.flamegpu.com/
   :Docs: http://docs.flamegpu.com/
   :GitHub: https://github.com/FLAMEGPU/FLAMEGPU

`FLAMEGPU <http://www.flamegpu.com>`_ is a multi agent simulation framework developed at the University of Sheffield.
It uses GPUs to accelerate the simulation of multiagent systems but abstracts the GPU architecture away from users so that they can write models using a high level of abstraction (without having to write GPU code).


Installing FLAMEGPU on Iceberg
------------------------------

To install FLAME GPU you should checkout the latest master branch of the code which has Linux support ::

    git clone https://github.com/FLAMEGPU/FLAMEGPU.git

To run the FLAME GPU examples you will need to be on a GPU node. You can start an interactive session using ::

    qrshx -l gpu=1 -l gpu_arch=nvidia-k40m -l rmem=13G

You will then need to load the relevant modules ::

    module load libs/cuda/8.0.44
    module load compilers/gcc/4.9.2


You can now navigate to the FLAME GPU examples folder and build the examples. e.g. to build all the examples in console mode ::

    cd FLAMEGPU/examples
    make console

Building the full suite of examples can take a while, instead you may wish to build individual examples. eg. to build a single example in console mode ::

    cd FLAMEGPU/examples/GameOfLife
    make console

You can now run the console versions of the example models by navigating to ``FLAMEGPU/bin/x64`` and calling the appropriate shell script. e.g. from the ``FLAMEGPU/examples/GameOfLife`` directory ::

    ../../bin/linux-x64/GameOfLife_console.sh

or by manually running the executable with the appropriate command line arguments. e.g. from the ``FLAMEGPU/examples/GameOfLife`` directory ::

   ../../bin/linux-x64/Release_Console/GameOfLife iterations/0.xml 1



FLAME GPU will run the example for one iteration and output the model state to a file ``1.xml`` in the model's iterations directory.

Visualisation currently does not work with X forwarding as FLAME GPU uses complex rendering techniques which are not supported.


For more information please see the `FLAME GPU Documentation <http://docs.flamegpu.com>`_. 
Additional information on compiling FLAME GPU examples can be found via :: 
    
    make help

And additional information on the command line arguments via :: 

    ../../bin/linux-x64/Release_Console/GameOfLife --help
