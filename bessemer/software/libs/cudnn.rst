.. _cudnn_bessemer:

cuDNN
=====

.. sidebar:: cuDNN

   :Dependencies: CUDA, gcc
   :URL: https://developer.nvidia.com/cudnn
   :Documentation: https://developer.nvidia.com/cudnn


The NVIDIA CUDA Deep Neural Network library (cuDNN) is
a GPU-accelerated library of primitives for deep neural networks.
cuDNN provides highly tuned implementations for standard routines such
as forward and backward convolution, pooling, normalization, and activation layers.
cuDNN is part of the NVIDIA Deep Learning SDK.

Usage
-----

**Only GPU-enabled nodes are able to run the library.**

.. See** :ref:`GPUComputing_bessemer` **for more information on how to request a GPU-enabled node for an interactive session or job submission.**

Load the appropriate cuDNN version
(and implicitly load a :ref:`specified CUDA version <cuda_bessemer>`)
with one of the following commands:

.. code-block:: bash

   module load cuDNN/7.4.2.24-gcccuda-2019a
   module load cuDNN/7.4.2.24-CUDA-10.0.130

Installation notes
------------------

This section is primarily for administrators of the system.

The cuDNN library is only available to download through the `developer portal <https://developer.nvidia.com/cudnn>`_.
Installation ``.tgz`` files must be located in ``/usr/local/media/eb-srcs/c/cuDNN/`` before cuDNN can be installed via EasyBuild.

cuDNN/7.4.2.24-gcccuda-2019a
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Installed using Easybuild-provided ``cuDNN-7.4.2.24-gcccuda-2019a`` easyconfig.

Testing: ran the ``mnistCUDNN`` example (see *Examples* above) with CUDA 10.1 (``foss-2019a`` toolchain) on a V100 GPU; results: ::

    /usr/local/packages/staging/eb/CUDA/10.1.105-GCC-8.2.0-2.31.1/bin/nvcc -ccbin g++ -I/usr/local/packages/staging/eb/CUDA/10.1.105-GCC-8.2.0-2.31.1/include -IFreeImage/include  -m64    -gencode arch=compute_30,code=sm_30 -gencode arch=compute_35,code=sm_35 -gencode arch=compute_50,code=sm_50 -gencode arch=compute_53,code=sm_53 -gencode arch=compute_53,code=compute_53 -o fp16_dev.o -c fp16_dev.cu
    g++ -I/usr/local/packages/staging/eb/CUDA/10.1.105-GCC-8.2.0-2.31.1/include -IFreeImage/include   -o fp16_emu.o -c fp16_emu.cpp
    g++ -I/usr/local/packages/staging/eb/CUDA/10.1.105-GCC-8.2.0-2.31.1/include -IFreeImage/include   -o mnistCUDNN.o -c mnistCUDNN.cpp
    /usr/local/packages/staging/eb/CUDA/10.1.105-GCC-8.2.0-2.31.1/bin/nvcc -ccbin g++   -m64      -gencode arch=compute_30,code=sm_30 -gencode arch=compute_35,code=sm_35 -gencode arch=compute_50,code=sm_50 -gencode arch=compute_53,code=sm_53 -gencode arch=compute_53,code=compute_53 -o mnistCUDNN fp16_dev.o fp16_emu.o mnistCUDNN.o -I/usr/local/packages/staging/eb/CUDA/10.1.105-GCC-8.2.0-2.31.1/include -IFreeImage/include  -LFreeImage/lib/linux/x86_64 -LFreeImage/lib/linux -lcudart -lcublas -lcudnn -lfreeimage -lstdc++ -lm
    cudnnGetVersion() : 7402 , CUDNN_VERSION from cudnn.h : 7402 (7.4.2)
    Host compiler version : GCC 8.2.0
    There are 1 CUDA capable devices on your machine :
    device 0 : sms 80  Capabilities 7.0, SmClock 1380.0 Mhz, MemSize (Mb) 32480, MemClock 877.0 Mhz, Ecc=1, boardGroupID=0
    Using device 0

    Testing single precision
    Loading image data/one_28x28.pgm
    Performing forward propagation ...
    Testing cudnnGetConvolutionForwardAlgorithm ...
    Fastest algorithm is Algo 0
    Testing cudnnFindConvolutionForwardAlgorithm ...
    ^^^^ CUDNN_STATUS_SUCCESS for Algo 0: 0.016384 time requiring 0 memory
    ^^^^ CUDNN_STATUS_SUCCESS for Algo 2: 0.047104 time requiring 57600 memory
    ^^^^ CUDNN_STATUS_SUCCESS for Algo 1: 0.081920 time requiring 3464 memory
    ^^^^ CUDNN_STATUS_SUCCESS for Algo 7: 0.084992 time requiring 2057744 memory
    ^^^^ CUDNN_STATUS_SUCCESS for Algo 5: 0.117760 time requiring 203008 memory
    Resulting weights from Softmax:
    0.0000000 0.9999399 0.0000000 0.0000000 0.0000561 0.0000000 0.0000012 0.0000017 0.0000010 0.0000000 
    Loading image data/three_28x28.pgm
    Performing forward propagation ...
    Resulting weights from Softmax:
    0.0000000 0.0000000 0.0000000 0.9999288 0.0000000 0.0000711 0.0000000 0.0000000 0.0000000 0.0000000 
    Loading image data/five_28x28.pgm
    Performing forward propagation ...
    Resulting weights from Softmax:
    0.0000000 0.0000008 0.0000000 0.0000002 0.0000000 0.9999820 0.0000154 0.0000000 0.0000012 0.0000006 

    Result of classification: 1 3 5

    Test passed!

    Testing half precision (math in single precision)
    Loading image data/one_28x28.pgm
    Performing forward propagation ...
    Testing cudnnGetConvolutionForwardAlgorithm ...
    Fastest algorithm is Algo 0
    Testing cudnnFindConvolutionForwardAlgorithm ...
    ^^^^ CUDNN_STATUS_SUCCESS for Algo 0: 0.016384 time requiring 0 memory
    ^^^^ CUDNN_STATUS_SUCCESS for Algo 2: 0.044032 time requiring 28800 memory
    ^^^^ CUDNN_STATUS_SUCCESS for Algo 1: 0.046080 time requiring 3464 memory
    ^^^^ CUDNN_STATUS_SUCCESS for Algo 7: 0.055296 time requiring 2057744 memory
    ^^^^ CUDNN_STATUS_SUCCESS for Algo 5: 0.068608 time requiring 203008 memory
    Resulting weights from Softmax:
    0.0000001 1.0000000 0.0000001 0.0000000 0.0000563 0.0000001 0.0000012 0.0000017 0.0000010 0.0000001 
    Loading image data/three_28x28.pgm
    Performing forward propagation ...
    Resulting weights from Softmax:
    0.0000000 0.0000000 0.0000000 1.0000000 0.0000000 0.0000714 0.0000000 0.0000000 0.0000000 0.0000000 
    Loading image data/five_28x28.pgm
    Performing forward propagation ...
    Resulting weights from Softmax:
    0.0000000 0.0000008 0.0000000 0.0000002 0.0000000 1.0000000 0.0000154 0.0000000 0.0000012 0.0000006 

    Result of classification: 1 3 5

    Test passed!

7.4.2.24-CUDA-10.0.130
^^^^^^^^^^^^^^^^^^^^^^

Installed using Easybuild-provided ``cuDNN-7.4.2.24-CUDA-10.0.130`` easyconfig.

Testing: ran the ``mnistCUDNN`` example (see *Examples* above) with CUDA 10.0 on a V100 GPU; results: ::

    /usr/local/packages/staging/eb/CUDA/10.0.130/bin/nvcc -ccbin g++ -I/usr/local/packages/staging/eb/CUDA/10.0.130/include -IFreeImage/include  -m64    -gencode arch=compute_30,code=sm_30 -gencode arch=compute_35,code=sm_35 -gencode arch=compute_50,code=sm_50 -gencode arch=compute_53,code=sm_53 -gencode arch=compute_53,code=compute_53 -o fp16_dev.o -c fp16_dev.cu
    g++ -I/usr/local/packages/staging/eb/CUDA/10.0.130/include -IFreeImage/include   -o fp16_emu.o -c fp16_emu.cpp
    g++ -I/usr/local/packages/staging/eb/CUDA/10.0.130/include -IFreeImage/include   -o mnistCUDNN.o -c mnistCUDNN.cpp
    /usr/local/packages/staging/eb/CUDA/10.0.130/bin/nvcc -ccbin g++   -m64      -gencode arch=compute_30,code=sm_30 -gencode arch=compute_35,code=sm_35 -gencode arch=compute_50,code=sm_50 -gencode arch=compute_53,code=sm_53 -gencode arch=compute_53,code=compute_53 -o mnistCUDNN fp16_dev.o fp16_emu.o mnistCUDNN.o -I/usr/local/packages/staging/eb/CUDA/10.0.130/include -IFreeImage/include  -LFreeImage/lib/linux/x86_64 -LFreeImage/lib/linux -lcudart -lcublas -lcudnn -lfreeimage -lstdc++ -lm
    cudnnGetVersion() : 7402 , CUDNN_VERSION from cudnn.h : 7402 (7.4.2)
    Host compiler version : GCC 4.8.5
    There are 1 CUDA capable devices on your machine :
    device 0 : sms 80  Capabilities 7.0, SmClock 1380.0 Mhz, MemSize (Mb) 32480, MemClock 877.0 Mhz, Ecc=1, boardGroupID=0
    Using device 0

    Testing single precision
    Loading image data/one_28x28.pgm
    Performing forward propagation ...
    Testing cudnnGetConvolutionForwardAlgorithm ...
    Fastest algorithm is Algo 0
    Testing cudnnFindConvolutionForwardAlgorithm ...
    ^^^^ CUDNN_STATUS_SUCCESS for Algo 0: 0.017408 time requiring 0 memory
    ^^^^ CUDNN_STATUS_SUCCESS for Algo 2: 0.048128 time requiring 57600 memory
    ^^^^ CUDNN_STATUS_SUCCESS for Algo 1: 0.083936 time requiring 3464 memory
    ^^^^ CUDNN_STATUS_SUCCESS for Algo 7: 0.088096 time requiring 2057744 memory
    ^^^^ CUDNN_STATUS_SUCCESS for Algo 4: 0.111616 time requiring 207360 memory
    Resulting weights from Softmax:
    0.0000000 0.9999399 0.0000000 0.0000000 0.0000561 0.0000000 0.0000012 0.0000017 0.0000010 0.0000000 
    Loading image data/three_28x28.pgm
    Performing forward propagation ...
    Resulting weights from Softmax:
    0.0000000 0.0000000 0.0000000 0.9999288 0.0000000 0.0000711 0.0000000 0.0000000 0.0000000 0.0000000 
    Loading image data/five_28x28.pgm
    Performing forward propagation ...
    Resulting weights from Softmax:
    0.0000000 0.0000008 0.0000000 0.0000002 0.0000000 0.9999820 0.0000154 0.0000000 0.0000012 0.0000006 

    Result of classification: 1 3 5

    Test passed!

    Testing half precision (math in single precision)
    Loading image data/one_28x28.pgm
    Performing forward propagation ...
    Testing cudnnGetConvolutionForwardAlgorithm ...
    Fastest algorithm is Algo 0
    Testing cudnnFindConvolutionForwardAlgorithm ...
    ^^^^ CUDNN_STATUS_SUCCESS for Algo 0: 0.014336 time requiring 0 memory
    ^^^^ CUDNN_STATUS_SUCCESS for Algo 7: 0.052224 time requiring 2057744 memory
    ^^^^ CUDNN_STATUS_SUCCESS for Algo 2: 0.052256 time requiring 28800 memory
    ^^^^ CUDNN_STATUS_SUCCESS for Algo 1: 0.056320 time requiring 3464 memory
    ^^^^ CUDNN_STATUS_SUCCESS for Algo 4: 0.067584 time requiring 207360 memory
    Resulting weights from Softmax:
    0.0000001 1.0000000 0.0000001 0.0000000 0.0000563 0.0000001 0.0000012 0.0000017 0.0000010 0.0000001 
    Loading image data/three_28x28.pgm
    Performing forward propagation ...
    Resulting weights from Softmax:
    0.0000000 0.0000000 0.0000000 1.0000000 0.0000000 0.0000720 0.0000000 0.0000000 0.0000000 0.0000000 
    Loading image data/five_28x28.pgm
    Performing forward propagation ...
    Resulting weights from Softmax:
    0.0000000 0.0000008 0.0000000 0.0000002 0.0000000 1.0000000 0.0000154 0.0000000 0.0000012 0.0000006 

    Result of classification: 1 3 5

    Test passed!

