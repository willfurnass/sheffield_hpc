.. _nvidia_compiler_sharc:

NVCC - Nvidia CUDA Compiler
===========================



To compile GPU code using the NVIDIA compiler, nvcc, first start an interactive session (see :ref:`GPUInteractive_sharc`).
Next, you need to set up the compiler environment via one of the following module statements: ::

    module load libs/cuda/8.0.44
    module load libs/cuda/7.5.18

depending on the version of CUDA you intend to use. This makes the ``nvcc`` CUDA compiler available.

**Important** To compile CUDA programs you also need a compatible version of the :ref:`gcc_sharc`.  As of version 8.0.44, CUDA is compatible with GCC versions:

* greater than or equal to 4.7.0 (to allow for the use of c++11 features) and
* less than 5.0.0

It is therefore recommended that you load the most recent 4.x version of GCC when building CUDA programs on ShARC: ::

    module load compilers/gcc/4.9.2

An example of the use of ``nvcc``::

    nvcc filename.cu

will compile the CUDA program contained in the file ``filename.cu``.


GPU Code Generation Options
---------------------------

To achieve the best possible performance whilst being portable, GPU code should be generated for the architecture(s) it will be executed upon.

This is controlled by specifying ``-gencode`` arguments to NVCC which, unlike the ``-arch`` and ``-code`` arguments, allows for fatbinary executables that are optimised for multiple device architectures.

Each ``-gencode`` argument requires two values, the *virtual architecture* and *real architecture*, for use in NVCC's `two-stage compilation <http://docs.nvidia.com/cuda/cuda-compiler-driver-nvcc/index.html#virtual-architectures>`_.
I.e. ``-gencode=arch=compute_60,code=sm_60`` specifies a virtual architecture of ``compute_60`` and real architecture ``sm_60``.

To support future hardware of higher compute capability, an additional ``-gencode`` argument can be used to enable Just in Time (JIT) compilation of embedded intermediate PTX code. This argument should use the highest virtual architecture specified in other gencode arguments for both the ``arch`` and ``code``, i.e. ``-gencode=arch=compute_60,code=sm_60``.

The minimum specified virtual architecture must be less than or equal to the `Compute Capability <https://developer.nvidia.com/cuda-gpus>`_ of the GPU used to execute the code.

Iceberg contains Telsa M2070 and Tesla K40m GPUs, which are compute capability 20 and 35 respectively.
To build a CUDA application which targets any GPU on Iceberg, use the following ``-gencode`` arguments: ::

    nvcc filename.cu
    -gencode=arch=compute_20,code=sm_20
    -gencode=arch=compute_35,code=sm_35
    -gencode=arch=compute_35,code=compute_35

ShARC contains Tesla K80 GPUs and Telsa P100 GPUs, which are compute capability 37 and 60 respectively.
To build a CUDA application which targets any GPU on ShARC, use the following ``-gencode`` arguments: ::

    nvcc filename.cu
    -gencode=arch=compute_37,code=sm_37
    -gencode=arch=compute_60,code=sm_60
    -gencode=arch=compute_60,code=compute_60

To build a CUDA application which targets any GPU on Icerberg or ShARC, use the following ``-gencode`` arguments: ::

    nvcc filename.cu
    -gencode=arch=compute_20,code=sm_20
    -gencode=arch=compute_35,code=sm_35
    -gencode=arch=compute_37,code=sm_37
    -gencode=arch=compute_60,code=sm_60
    -gencode=arch=compute_60,code=compute_60
    

Further details of these compiler flags can be found in the `NVCC Documentation <http://docs.nvidia.com/cuda/cuda-compiler-driver-nvcc/index.html#options-for-steering-gpu-code-generation>`_, 
along with details of the supported `virtual architectures <http://docs.nvidia.com/cuda/cuda-compiler-driver-nvcc/index.html#virtual-architecture-feature-list>`_ and `real architectures <http://docs.nvidia.com/cuda/cuda-compiler-driver-nvcc/index.html#gpu-feature-list>`_.

.. note:: SM 20 and SM 21 are deprecated in CUDA 8.0.

  If you attempt to build SM 20 or SM 21 code using CUDA 8.0, a warning will be raised at compile time. 

.. warning:: SM 20 and SM 21 are removed in CUDA 8.0.

  It is **not possible** to build SM 20 or SM 21 code using CUDA 9.0 or above. 
  
  If you are using CUDA 9.0 or greater remove ``-gencode`` arguments containing ``compute`` or ``sm`` values of 20 and 21 from the above examples.
  
  It is also not possible to execute code generated using CUDA 9.0 or above on the Tesla M2070s in Iceberg. 

