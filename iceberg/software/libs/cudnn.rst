.. _cudnn:

cuDNN
=====

.. sidebar:: cudNN

  
   :Dependencies: CUDA, gcc
   :URL: https://developer.nvidia.com/cudnn
   :Documentation: https://developer.nvidia.com/cudnn


The NVIDIA CUDA Deep Neural Network library (cuDNN) is a GPU-accelerated library of primitives for deep neural networks. cuDNN provides highly tuned implementations for standard routines such as forward and backward convolution, pooling, normalization, and activation layers. cuDNN is part of the NVIDIA Deep Learning SDK.

Usage
-----

Currently cuDNN 5.0 is available for CUDA versions 8.0.x and 7.5.x. An appropriate **CUDA module is loaded automatically** so there's no need for a separate CUDA module load call.

Load the appropriate cuDNN version with one of the following commands: ::

    module load libs/binlibs/cudnn/cuda8.0/cudnn5.0
    module load libs/binlibs/cudnn/cuda7.5/cudnn5.0    


Installation notes
------------------

This section is primarily for administrators of the system.

- The cuDNN library is only available to download through the `developer portal <https://developer.nvidia.com/cudnn>`_.
- Install script: `install_cudnn5.0_for_cuda7.5_cuda8.0.sh <https://github.com/rcgsheffield/sheffield_hpc/tree/master/iceberg/software/install_scripts/libs/binlibs/cudnn/install_cudnn5.0_for_cuda7.5_cuda8.0.sh>`_
- Module file: 
	- `cuDNN 5.0 for CUDA 7.5 <https://github.com/rcgsheffield/sheffield_hpc/tree/master/iceberg/software/modulefiles/libs/binlibs/cudnn/cuda7.5/cudnn5.0>`_
	- `cuDNN 5.0 for CUDA 8.0 <https://github.com/rcgsheffield/sheffield_hpc/tree/master/iceberg/software/modulefiles/libs/binlibs/cudnn/cuda8.0/cudnn5.0>`_





