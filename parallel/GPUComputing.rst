.. _GPUIntro:

GPU Computing
=============

Graphics Processing Units (GPUs) were, as the name suggests, originally designed for the efficient processing of graphics. Over time, they were developed into systems that were capable of performing general purpose computing which is why some people refer to modern GPUs as GP-GPUs (General Purpose Graphical Processing Units).

Graphics processing tpically involves relatively simple computations that need to applied to millions of on-screen pixels in parallel. As such, GPUs tend to be very quick and efficient at computing certain types of parallel workloads.

The `GPUComputing@Sheffield website
<http://gpucomputing.shef.ac.uk/>`_ aims to facilitate the use of GPU computing within University of Sheffield research by providing resources for training and purchasing of equipment as well as providing a network of GPU users and research projects within the University.

.. _GPUCommunity:

GPU Community and NVIDIA Research Centre
----------------------------------------
The University of Sheffield has been officially affiliated with `NVIDIA
<https://research.nvidia.com/>`_ since 2011 as an `NVIDIA CUDA Research Centre
<https://developer.nvidia.com/academia/centers/university-sheffield-cuda-research-center>`_. As such NVIDIA offer us some benefits as a research institution including discounts on hardware, technical liaisons, online training and priority seed hardware for new GPU architectures. For first access to hardware, training and details of upcoming events, discussions and help please join the `GPUComputing google group
<https://groups.google.com/a/sheffield.ac.uk/forum/#!forum/gpucomputing>`_.

.. _GPUAccess:

Requesting access to GPU facilities
-----------------------------------

In order to ensure that the nodes hosting the GPUs run only GPU related tasks, we have defined a special project-group for accessing these nodes. If you wish to take advantage of the GPU processors, please contact research-it@sheffield.ac.uk asking to join the GPU project group.

Any Iceberg/ShARC user can apply to join this group. However, because our GPU resources are limited we will need to discuss the needs of the user and obtain consent from the project leader before allowing usage.

.. _GPUInteractive:

Interactive use of the GPUs
---------------------------

Once you are included in the GPU project group you may start using the GPU enabled nodes interactively by typing: ::

        qsh -l gpu=1 -P gpu

the ``-l gpu=`` parameter determines how many GPUs you are requesting. Currently, the maximum number of GPUs allowed per job is set to 4, i.e. you cannot exceed ``-l gpu=4``. Most jobs will only make use of one GPU.

**On Iceberg**, if your job requires selecting the type of GPU hardware, one of the following two optional parameters can be used to make that choice ::

	qsh -l gpu_arch=nvidia-m2070 -P gpu
	qsh -l gpu_arch=nvidia-k40m -P gpu

Interactive sessions provide you with 2 Gigabytes of CPU RAM by default which is significantly less than the amount of GPU RAM available. This can lead to issues where your session has insufficient CPU RAM to transfer data to and from the GPU. As such, it is recommended that you request enough CPU memory to communicate properly with the GPU ::

  qsh -l gpu_arch=nvidia-m2070 -P gpu -l rmem=7G
  qsh -l gpu_arch=nvidia-k40m -P gpu -l rmem=13G

The above will give you 1GB more CPU RAM than GPU RAM for each of the respective GPU architectures.

**On ShARC**, Nvidia K80 cards and P100s (in DGX-1) are available but access is provided using a queue based system instead so there is no need to call ``-l gpu_arch=``. Instead, use  ``-l gpu=`` together with ``-q your_queue_name.q``, e.g. : ::

  qsh -l gpu=1 -P your_project_name -q your_queue_name.q


.. _GPUJobs:

Submitting batch GPU jobs
-------------------------

To run batch jobs on gpu nodes, edit your jobfile to include a request for GPUs, e.g. for a single GPU ::

  #!/bin/bash
  #$ -l gpu=1 -P gpu


You can also use the the ``gpu_arch`` discussed aboved to target a specific GPU model ::

  #!/bin/bash
  #$ -l gpu_arch=nvidia-m2070 -P gpu


.. _GPUResources_iceberg:

Iceberg GPU Resources
---------------------

Hardware
^^^^^^^^

Iceberg currently contains 16 GPU units:

* 8 Nvidia Tesla Kepler K40M GPU units. Each unit contains 2880 CUDA cores, 12GB of memory and is capable of up to 1.43 Teraflops of double precision compute power.
* 8 Nvidia Tesla Fermi M2070 GPU units. Each unit contains 448 CUDA cores, 6GB of memory and is capable of up to 515 Gigaflops of double precision compute power.

GPU-enabled Software
^^^^^^^^^^^^^^^^^^^^

* Applications
    * :ref:`ansys_iceberg`
    * :ref:`maple_iceberg`
    * :ref:`matlab_iceberg`
    * :ref:`theano_iceberg`
* Libraries
    * :ref:`cuda_iceberg`
    * :ref:`cudnn_iceberg`
* Development Tools
    * :ref:`PGI Compilers`
    * :ref:`nvidia_compiler_iceberg`


.. _GPUResources_sharc:

ShARC GPU Resources
-------------------

Hardware
^^^^^^^^

ShARC currently contains 8 GPU units:

* 8 Nvidia Tesla Kepler K80 GPU units. Each unit contains 4992 CUDA cores, 24GB of memory and is capable of up to 2.91 Teraflops of double precision compute power.

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
