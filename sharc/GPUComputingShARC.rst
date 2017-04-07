.. _GPUComputing_sharc:

Using GPUs on ShARC
===================

Requesting access to GPU facilities
-----------------------------------

**Public GPU nodes have now been made available to Iceberg and ShARC users, these can be be used without acquiring extra permission.**

Research groups also have an option to purchase and add nodes to the ShARC cluster to be managed by CiCs (contact research-it@sheffield.ac.uk for more information). For these nodes (e.g. :ref:`dgx1_com_groupnodes_sharc`), permission from the group leader is required for access.

.. _GPUInteractive_sharc:

Interactive use of the GPUs
---------------------------

Once you are included in the GPU project group you may start using the GPU enabled nodes interactively by typing: ::

        qsh -l gpu=1

the ``-l gpu=`` parameter determines how many GPUs you are requesting. Currently, the maximum number of GPUs allowed per job is set to 4, i.e. you cannot exceed ``-l gpu=4``. Most jobs will only make use of one GPU.

Interactive sessions provide you with 2 Gigabytes of CPU RAM by default which is significantly less than the amount of GPU RAM available. This can lead to issues where your session has insufficient CPU RAM to transfer data to and from the GPU. As such, it is recommended that you request enough CPU memory to communicate properly with the GPU ::

  #Nvidia K80 GPU has 24GB of RAM
  qsh -l gpu=1 -l rmem=25G

The above will give you 1GB more CPU RAM than the 24GB of GPU RAM available on the Nvidia K80.


.. _GPUJobs_sharc:

Submitting batch GPU jobs
-------------------------

To run batch jobs on gpu nodes, edit your jobfile to include a request for GPUs, e.g. for a single GPU ::

  #!/bin/bash
  #$ -l gpu=1



.. _GPUResources_sharc:

ShARC GPU Resources
-------------------

Hardware
^^^^^^^^

**ShARC currently contains 8 publicly available GPU units:**

* 8 Nvidia Tesla Kepler K80 GPU units. Each unit contains 4992 CUDA cores, 24GB of memory and is capable of up to 2.91 Teraflops of double precision compute power.

**Research group-specific GPU nodes:**

* :ref:`dgx1_com_groupnodes_sharc`


GPU-enabled Software
^^^^^^^^^^^^^^^^^^^^

* Applications
    * :ref:`caffe_sharc`
    * :ref:`matlab_sharc`
    * :ref:`theano_sharc`
    * :ref:`tensorflow_sharc`
    * :ref:`torch_sharc`
* Libraries
    * :ref:`cuda_sharc`
    * :ref:`cudnn_sharc`
* Development Tools
    * :ref:`PGI Compilers_sharc`
    * :ref:`nvidia_compiler_sharc`
