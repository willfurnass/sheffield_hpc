.. _cudnn_sharc:

cuDNN
=====

.. sidebar:: cuDNN


   :Dependencies: CUDA, gcc
   :URL: https://developer.nvidia.com/cudnn
   :Documentation: https://developer.nvidia.com/cudnn


The NVIDIA CUDA Deep Neural Network library (cuDNN) is a GPU-accelerated library of primitives for deep neural networks. cuDNN provides highly tuned implementations for standard routines such as forward and backward convolution, pooling, normalization, and activation layers. cuDNN is part of the NVIDIA Deep Learning SDK.

Usage
-----

**Only GPU-enabled nodes are able to run the library. See** :ref:`GPUComputing_sharc` **for more information on how to request a GPU-enabled node for an interactive session or job submission.**

Load the appropriate cuDNN version (**and implicitly load a specific CUDA version**) with one of the following commands: ::

    module load libs/cudnn/7.3.1.20/binary-cuda-9.0.176
    module load libs/cudnn/7.0/binary-cuda-9.1.85
    module load libs/cudnn/7.0/binary-cuda-8.0.44
    module load libs/cudnn/6.0/binary-cuda-8.0.44
    module load libs/cudnn/5.1/binary-cuda-8.0.44
    module load libs/cudnn/5.1/binary-cuda-7.5.18

Examples
--------

Examples are provided with ``libs/cudnn/7.3.1.20/binary-cuda-9.0.176``: 

.. code_block:: bash

   # start an interactive session on a GPU node
   qrshx -l gpu=1  
   # Copy the cuDNN 7 examples to a temporary directory
   cp -r $CUDNN_HOME/src/cudnn_samples_v7 $TMPDIR
   # Compile an example
   cd mnistCUDNN
   make
   # Run the example
   ./mnistCUDNN

Installation notes
------------------

This section is primarily for administrators of the system.

The cuDNN library is only available to download through the `developer portal <https://developer.nvidia.com/cudnn>`_.  Installation ``.tgz`` files are located in ``/usr/local/media/protected/cudnn``.

Version 7.3.1.20
^^^^^^^^^^^^^^^^

- Install script: :download:`install_cudnn7.3.1.20_for_cuda9.0.sh </sharc/software/install_scripts/libs/cudnn/install_7.0_for_cuda9.0.sh>`
- :download:`Module file for CUDA 9.0 </sharc/software/modulefiles/libs/cudnn/7.3.1.20/binary-cuda-9.0.176>`
- Testing: ran the ``mnistCUDNN`` example (see *Examples* above) with CUDA 9.0 on a V100 GPU; results: ::

  [te1st@sharc-node168 mnistCUDNN]$ ./mnistCUDNN 
  cudnnGetVersion() : 7301 , CUDNN_VERSION from cudnn.h : 7301 (7.3.1)
  Host compiler version : GCC 4.8.5
  There are 2 CUDA capable devices on your machine :
  device 0 : sms 80  Capabilities 7.0, SmClock 1380.0 Mhz, MemSize (Mb) 16160, MemClock 877.0 Mhz, Ecc=1, boardGroupID=0
  device 1 : sms 80  Capabilities 7.0, SmClock 1380.0 Mhz, MemSize (Mb) 16160, MemClock 877.0 Mhz, Ecc=1, boardGroupID=1
  Using device 0
  
  Testing single precision
  Loading image data/one_28x28.pgm
  Performing forward propagation ...
  Testing cudnnGetConvolutionForwardAlgorithm ...
  Fastest algorithm is Algo 0
  Testing cudnnFindConvolutionForwardAlgorithm ...
  ^^^^ CUDNN_STATUS_SUCCESS for Algo 0: 0.109600 time requiring 0 memory
  ^^^^ CUDNN_STATUS_SUCCESS for Algo 2: 0.161792 time requiring 57600 memory
  ^^^^ CUDNN_STATUS_SUCCESS for Algo 1: 0.183296 time requiring 3464 memory
  ^^^^ CUDNN_STATUS_SUCCESS for Algo 5: 0.198656 time requiring 203008 memory
  ^^^^ CUDNN_STATUS_SUCCESS for Algo 7: 0.200704 time requiring 2057744 memory
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
  ^^^^ CUDNN_STATUS_SUCCESS for Algo 0: 0.048128 time requiring 0 memory
  ^^^^ CUDNN_STATUS_SUCCESS for Algo 1: 0.089088 time requiring 3464 memory
  ^^^^ CUDNN_STATUS_SUCCESS for Algo 7: 0.097280 time requiring 2057744 memory
  ^^^^ CUDNN_STATUS_SUCCESS for Algo 2: 0.098272 time requiring 28800 memory
  ^^^^ CUDNN_STATUS_SUCCESS for Algo 4: 0.132096 time requiring 207360 memory
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

Version 7.0
^^^^^^^^^^^

- Install script: :download:`install_cudnn7.0_for_cuda8.0_cuda9.1.sh </sharc/software/install_scripts/libs/cudnn/install_7.0_for_cuda_8.0_cuda_9.1.sh>`
- :download:`Module file for CUDA 9.1 </sharc/software/modulefiles/libs/cudnn/7.0/binary-cuda-9.1.85>`
- :download:`Module file for CUDA 8.0 </sharc/software/modulefiles/libs/cudnn/7.0/binary-cuda-8.0.44>`

Version 6.0
^^^^^^^^^^^

- Install script: :download:`install_cudnn6.0_for_cuda8.0.sh </sharc/software/install_scripts/libs/cudnn/install_6.0_for_cuda_8.0.sh>`
- :download:`Module file </sharc/software/modulefiles/libs/cudnn/6.0/binary-cuda-8.0.44>`

Version 5.1
^^^^^^^^^^^

- Install script: :download:`install_cudnn5.1_for_cuda7.5_cuda8.0.sh </sharc/software/install_scripts/libs/cudnn/install_5.1_for_cuda_7.5_cuda_8.0.sh>`
- :download:`Module file for CUDA 8.0 </sharc/software/modulefiles/libs/cudnn/5.1/binary-cuda-8.0.44>`
- :download:`Module file for CUDA 7.5 </sharc/software/modulefiles/libs/cudnn/5.1/binary-cuda-7.5.18>`





