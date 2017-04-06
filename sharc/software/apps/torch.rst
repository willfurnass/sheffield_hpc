.. _torch_sharc:

Torch
=====

.. sidebar:: Torch

   :URL: http://torch.ch/

Torch is a scientific computing framework with wide support for machine learning algorithms that puts GPUs first. It is easy to use and efficient, thanks to an easy and fast scripting language, LuaJIT, and an underlying C/CUDA implementation.

**Additional permissions are needed to use GPUs on Iceberg/ShARC. See** :ref:`GPUComputing_sharc` **for more information.**

Torch Singularity Images
------------------------

Singularity images are self-contained virtual machines similar to Docker. For more information on Singularity and how to use the images, see :ref:`singularity_sharc`.

The following Singularity images are available on ShARC and can also be downloaded for use on your local machine:

* GPU Torch 7, Ubuntu 16.04, CUDA 8, cuDNN 5.0, GCC 5.4.0
    * Image path: ``/usr/local/packages/singularity/torch/v7-GPU-Ubuntu16.04-CUDA8-cudNN5.0.img``
    * Def file: `Caffe GPU </sharc/software/apps/singularity/torch_gpu.def>`

To get a bash terminal in to an image for example, use the command: ::

  singularity exec /usr/local/packages/singularity/torch/v7-GPU-Ubuntu16.04-CUDA8-cudNN5.0.img /bin/bash

The ``exec`` command can also be used to call any command/script inside the image e.g. ::

  singularity exec /usr/local/packages/singularity/torch/v7-GPU-Ubuntu16.04-CUDA8-cudNN5.0.img th

The paths ``/fastdata``, ``/data``, ``/home``, ``/scratch``, ``/shared`` are automatically mounted to your ShARC filestore directories. For GPU-enabled images the ``/nvlib`` and ``/nvbin`` is mounted to the correct Nvidia driver version for the node that you're using.


Using the Torch Module
----------------------

The following is an instruction on how to use the Torch module.

First request an interactive session, e.g. with :ref:`qrshx`. To use GPUs see :ref:`GPUInteractive_sharc`.

Load the Torch module which also loads anaconda 3.4, CUDA 8.0, cuDNN 5.1 and GCC 4.9.4. ::

	module load apps/torch/nvidia-7/gcc-4.9.4-cuda-8.0-cudnn-5.1-conda-3.4-TESTING

On the DGX-1, load the NCCL library optimised for the hardware: ::

	moudle libs/nccl/dgx-1/binary-cuda-8.0

On any other node use a generic NCCL build: ::

	moudle load libs/nccl/generic/gcc-4.9.4-cuda-8.0


(Optional if you plan to interface with python) Create a conda environment to load relevant modules on your local user account and activate it ::

	conda create -n torch python=3.5
	source activate torch



Every Session Afterwards and in Your Job Scripts
------------------------------------------------

In the interactive session or your batch script, load the relevant modules and (optionally) activate your conda environment ::

	module load apps/torch/nvidia-7/gcc-4.9.4-cuda-8.0-cudnn-5.1-conda-3.4-TESTING

	#Optional
	source activate torch
