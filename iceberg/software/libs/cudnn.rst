.. _iceberg_cudnn:

cuDNN
=====

.. sidebar:: cudNN

  
   :Dependencies: CUDA, gcc
   :URL: https://developer.nvidia.com/cudnn
   :Documentation: https://developer.nvidia.com/cudnn


The NVIDIA CUDA Deep Neural Network library (cuDNN) is a GPU-accelerated library of primitives for deep neural networks. cuDNN provides highly tuned implementations for standard routines such as forward and backward convolution, pooling, normalization, and activation layers. cuDNN is part of the NVIDIA Deep Learning SDK.

Usage
-----

Currently cuDNN 5.1 is available for CUDA versions 8.0.44 and 7.5.18. An appropriate **CUDA module is loaded automatically** so there's no need for a separate CUDA module load call.

Load the appropriate cuDNN version with one of the following commands: ::

    module load libs/binlibs/cudnn/5.1-cuda-8.0.44
    module load libs/binlibs/cudnn/5.1-cuda-7.5.18    


Installation notes
------------------

This section is primarily for administrators of the system.

- The cuDNN library is only available to download through the `developer portal <https://developer.nvidia.com/cudnn>`_.
- Installation
	- Install script: :download:`install_cudnn5.1_for_cuda7.5_cuda8.0.sh </iceberg/software/install_scripts/libs/binlibs/cudnn/install_cudnn5.1_for_cuda7.5_cuda8.0.sh>`
	- Installation ``.tgz`` files are located in ``/usr/local/media/protected/cudnn``
- Module file: 
	- :download:`cuDNN 5.1 for CUDA 7.5 </iceberg/software/modulefiles/libs/binlibs/cudnn/5.1-cuda-7.5.18>`
	- :download:`cuDNN 5.1 for CUDA 8.0 </iceberg/software/modulefiles/libs/binlibs/cudnn/5.1-cuda-8.0.44>`





