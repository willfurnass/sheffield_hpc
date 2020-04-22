.. _`cuda_iceberg`:

CUDA
====

CUDA (*Compute Unified Device Architecture*) 
is a parallel computing platform and application programming interface (API) model created by NVIDIA.
It allows software developers to use a CUDA-enabled graphics processing unit (GPU) for general purpose processing, 
an approach known as *General Purpose GPU* (GPGPU) computing.

Usage
-----

You need to first request one or more GPUs within an :ref:`interactive session or batch job on a worker node <submit-queue>`.  
For example, to request a single GPU for an interactive session on a worker node:

.. code-block:: sh

   qrshx -l gpu=1

.. note:: See :ref:`GPUComputing_iceberg` for more information on how to request a GPU-enabled node for an interactive session or job submission. 

You then need to load a version of the CUDA library (and compiler).
There are several versions of the CUDA library available. 
As with much software installed on the cluster, 
versions of CUDA are activated via the :ref:`'module load' command<env_modules>`:

.. code-block:: sh

   module load libs/cuda/8.0.44
   module load libs/cuda/7.5.18
   module load libs/cuda/6.5.14
   module load libs/cuda/4.0.17
   module load libs/cuda/3.2.16

To then confirm which version of CUDA you are using:

.. code-block:: console

   $ nvcc --version
   nvcc: NVIDIA (R) Cuda compiler driver
   Copyright (c) 2005-2016 NVIDIA Corporation
   Built on Sun_Sep__4_22:14:01_CDT_2016
   Cuda compilation tools, release 8.0, V8.0.44

**Important** To compile CUDA programs you also need a compatible version of the :ref:`GCC compiler <gcc_iceberg>`:

* CUDA 7.x and 8.x: GCC >= 4.7.0 (to allow for the use of c++11 features) and < 5.0.0

Compiling a simple CUDA program
-------------------------------

An example of the use of ``nvcc`` (the CUDA compiler)

.. code-block:: sh

   nvcc filename.cu

will compile the CUDA program contained in the file ``filename.cu``.

Compiling the sample programs
-----------------------------

You do not need to be using a GPU-enabled node to compile the sample programs but you do need a GPU to run them.

In a ``qrshx`` session:

.. code-block:: sh

   # Load modules
   module load libs/cuda/8.0.44
   module load compilers/gcc/4.9.2

   # Copy CUDA samples to a local directory
   # It will create a directory called NVIDIA_CUDA-8.0_Samples/
   mkdir cuda_samples
   cd cuda_samples
   cp -r $CUDA_SDK .

   # Compile (this will take a while)
   cd NVIDIA_CUDA-8.0_Samples/
   make

The ``make`` command then runs the ``nvcc`` CUDA compiler and
generates a binary executable that you can then run on a node with
an NVIDIA GPU installed.

A basic test is to run one of the resulting binaries, ``deviceQuery``.

GPU Code Generation Options
---------------------------

To achieve the best possible performance whilst being portable, 
GPU code should be generated for the architecture(s) it will be executed upon.

This is controlled by specifying ``-gencode`` arguments to NVCC which, 
unlike the ``-arch`` and ``-code`` arguments, 
allows for 'fatbinary' executables that are optimised for multiple device architectures.

Each ``-gencode`` argument requires two values, 
the *virtual architecture* and *real architecture*, 
for use in NVCC's `two-stage compilation <https://docs.nvidia.com/cuda/cuda-compiler-driver-nvcc/index.html#virtual-architectures>`_.
I.e. ``-gencode=arch=compute_20,code=sm_20`` specifies a virtual architecture of ``compute_20`` and real architecture ``sm_20``.

To support future hardware of higher compute capability, 
an additional ``-gencode`` argument can be used to enable Just in Time (JIT) compilation of embedded intermediate PTX code. 
This argument should use the highest virtual architecture specified in other gencode arguments 
for both the ``arch`` and ``code``
i.e ``-gencode=arch=compute_20,code=compute_20``.

The minimum specified virtual architecture must be less than or equal to the `Compute Capability <https://developer.nvidia.com/cuda-gpus>`_ of the GPU used to execute the code.

Iceberg contains Telsa M2070 and Tesla K40m GPUs, 
which are compute capability 20 and 35 respectively.
To build a CUDA application which targets any GPU on Iceberg, 
use the following ``-gencode`` arguments:

.. code-block:: sh

   nvcc filename.cu \
      -gencode=arch=compute_20,code=sm_20 \
      -gencode=arch=compute_35,code=sm_35 \
      -gencode=arch=compute_35,code=compute_35

To build a CUDA application that runs on both Iceberg and ShARC see :ref:`cuda_sharc`.

Further details of these compiler flags can be found in the `NVCC Documentation <https://docs.nvidia.com/cuda/cuda-compiler-driver-nvcc/index.html#options-for-steering-gpu-code-generation>`_, 
along with details of the supported `virtual architectures <https://docs.nvidia.com/cuda/cuda-compiler-driver-nvcc/index.html#virtual-architecture-feature-list>`_ and `real architectures <https://docs.nvidia.com/cuda/cuda-compiler-driver-nvcc/index.html#gpu-feature-list>`_.

.. note:: SM 20 and SM 21 are deprecated in CUDA 8.0.

  If you attempt to build SM 20 or SM 21 code using CUDA 8.0, a warning will be raised at compile time. 

.. warning:: SM 35 is not available in CUDA 3.2.16 or CUDA 4.0.17

  If you wish to target the Tesla K40m GPUs please use CUDA 6.5.14 or later.

Documentation
-------------

* `CUDA Toolkit Documentation <https://docs.nvidia.com/cuda/index.html#axzz3uLoSltnh>`_
* `The power of C++11 in CUDA 7 <http://devblogs.nvidia.com/parallelforall/power-cpp11-cuda-7/>`_

Profiling using nvprof
----------------------

Note that ``nvprof``, NVIDIA's CUDA profiler, 
cannot write output to the ``/fastdata`` filesystem.

This is because the profiler's output is a `SQLite <https://www.sqlite.org/>`__ database 
and SQLite requires a filesystem that supports file locking
but file locking is not enabled on the (`Lustre <http://lustre.org/>`__) filesystem mounted on ``/fastdata`` 
(for performance reasons). 

CUDA Training
-------------

`GPUComputing@sheffield <http://gpucomputing.shef.ac.uk>`_ provides 
a self-paced `introduction to CUDA <http://gpucomputing.shef.ac.uk/education/cuda/>`_ training course.

Determining the NVIDIA Driver version
-------------------------------------

Run the command:

.. code-block:: sh

   cat /proc/driver/nvidia/version

Example output is: ::

   NVRM version: NVIDIA UNIX x86_64 Kernel Module  384.81  Wed Aug 17 22:24:07 PDT 2016
   GCC version:  gcc version 4.4.7 20120313 (Red Hat 4.4.7-17) (GCC)

Installation notes
------------------

These are primarily for system administrators.

Device driver
^^^^^^^^^^^^^

The NVIDIA device driver is installed and configured using the ``/etc/init.d/uos-nvidia`` service.

This service does the following at boot time:

- Check the device driver version and uninstall it then reinstall the target version if required;
- Load the ``nvidia`` kernel module;
- Create several *device nodes* in ``/dev/``.

The NVIDIA device driver is currently version 384.81.  The driver installer provides OpenGL libraries.

CUDA 8.0.44
^^^^^^^^^^^

#. The CUDA toolkit binaries and samples were installed using a binary ``.run`` file:

   .. code-block:: sh

      cuda_vers="8.0.44"
      prefix="/usr/local/packages/libs/CUDA/binlibs${cuda_vers}"
      mkdir -m 2775 -p $prefix
      chown ${USER}:app-admins $prefix
      cd /usr/local/media/nvidia/
      chmod +x cuda_${cuda_vers}_linux.run
      ./cuda_${cuda_vers}_linux.run --toolkit --toolkitpath=${prefix}/cuda \
                                    --samples --samplespath=${prefix}/samples \
                                    --no-opengl-libs -silent

#. :download:`This modulefile </iceberg/software/modulefiles/libs/binlibs/cuda/8.0.44>` was installed as ``/usr/local/modulefiles/libs/cuda/8.0.44``

CUDA 7.5.18
^^^^^^^^^^^
**CUDA 7.5.18**

#. The CUDA toolkit binaries and samples were installed using a binary ``.run`` file as per version 8.0.44.
#. :download:`This modulefile </iceberg/software/modulefiles/libs/binlibs/cuda/7.5.18>` was installed as ``/usr/local/modulefiles/libs/cuda/7.5.18``

**Previous versions**

No install notes are available.
