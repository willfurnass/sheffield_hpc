.. _torch_sharc:

Torch
=====

.. sidebar:: Torch

   :URL: http://torch.ch/

Torch is a scientific computing framework with wide support for machine learning algorithms that puts GPUs first. It is easy to use and efficient, thanks to an easy and fast scripting language, LuaJIT, and an underlying C/CUDA implementation.

**Additional permissions are needed to use GPUs on Iceberg/ShARC. See** :ref:`GPUComputing_sharc` **for more information.**

Installation
------------

The following is an instruction on how to setup Torch on your local user account.

First request an interactive session, e.g. with :ref:`qrshx`. To use GPUs see :ref:`GPUInteractive_sharc`.

Load the Torch module which also loads anaconda 3.4, CUDA 8.0, cuDNN 5.1 and GCC 4.9.4. ::

	module load apps/torch/nvidia-7/gcc-4.9.4-cuda-8.0-cudnn-5.1-conda-3.4-TESTING

On the DGX-1, load the NCCL library optimised for the hardware: ::

	moudle libs/nccl/dgx-1/binary-cuda-8.0

On any other node use a generic NCCL build: ::

	moudle load libs/nccl/generic/gcc-4.9.4-cuda-8.0


(Optional if you plan to interace with python) Create a conda environment to load relevant modules on your local user account and activate it ::

	conda create -n torch python=3.5
	source activate torch



Every Session Afterwards and in Your Job Scripts
------------------------------------------------

In the interactive session or your batch script, load the relevant modules and (optionally) activate your conda environment ::

	module load apps/torch/nvidia-7/gcc-4.9.4-cuda-8.0-cudnn-5.1-conda-3.4-TESTING

	#Optional
	source activate torch
