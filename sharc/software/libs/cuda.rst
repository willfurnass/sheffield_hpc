.. _cuda_sharc:

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

.. note:: See :ref:`GPUComputing_sharc` for more information on how to request a GPU-enabled node for an interactive session or job submission. 

You then need to load a version of the CUDA library (and compiler).
There are several versions of the CUDA library available. 
As with much software installed on the cluster, 
versions of CUDA are activated via the :ref:`'module load' command<env_modules>`:

.. code-block:: sh

   module load libs/CUDA/10.2.89/binary
   module load libs/CUDA/10.1.243/binary
   module load libs/CUDA/10.0.130/binary
   module load libs/CUDA/9.1.85/binary
   module load libs/CUDA/9.0.176/binary
   module load libs/CUDA/8.0.44/binary
   module load libs/CUDA/7.5.18/binary

To then confirm which version of CUDA you are using:

.. code-block:: console

    $ nvcc --version
    nvcc: NVIDIA (R) Cuda compiler driver
    Copyright (c) 2005-2019 NVIDIA Corporation
    Built on Wed_Oct_23_19:24:38_PDT_2019
    Cuda compilation tools, release 10.2, V10.2.89

**Important** To compile CUDA programs you also need a compatible version of the :ref:`GCC compiler <gcc_sharc>`:

* CUDA 7.x and 8.x: GCC >= 4.7.0 (to allow for the use of c++11 features) and < 5.0.0
* CUDA 9.x: GCC < 7.0.0

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
   module load libs/CUDA/9.1.85/binary
   module load dev/gcc/4.9.4

   # Copy CUDA samples to a local directory
   # It will create a directory called NVIDIA_CUDA-9.1_Samples/
   mkdir cuda_samples
   cd cuda_samples
   cp -r $CUDA_SDK .

   # Compile (this will take a while)
   cd NVIDIA_CUDA-9.1_Samples/
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
I.e. ``-gencode=arch=compute_60,code=sm_60`` specifies a virtual architecture of ``compute_60`` and real architecture ``sm_60``.

To support future hardware of higher compute capability, 
an additional ``-gencode`` argument can be used to enable Just in Time (JIT) compilation of embedded intermediate PTX code. 
This argument should use the highest virtual architecture specified in other gencode arguments 
for both the ``arch`` and ``code``
i.e. ``-gencode=arch=compute_60,code=compute_60``.

The minimum specified virtual architecture must be less than or equal to the `Compute Capability <https://developer.nvidia.com/cuda-gpus>`_ of the GPU used to execute the code.

Public GPU nodes in ShARC contain Tesla K80 GPUs, which are compute capability 37.
To build a CUDA application which targets the public GPUS nodes, use the following ``-gencode`` arguments: 

.. code-block:: sh

   nvcc filename.cu \
      -gencode=arch=compute_37,code=sm_37 \
      -gencode=arch=compute_37,code=compute_37

ShARC also contains Tesla P100 GPUs and Tesla V100 GPUs in private nodes,
which are compute capability 60 and 70 respectively.
To build a CUDA application which targets any GPU on ShARC (either public or private), 
use the following ``-gencode`` arguments (for CUDA 9.0 and above):

.. code-block:: sh

   nvcc filename.cu \
      -gencode=arch=compute_37,code=sm_37 \
      -gencode=arch=compute_60,code=sm_60 \
      -gencode=arch=compute_70,code=sm_70 \
      -gencode=arch=compute_70,code=compute_70


Further details of these compiler flags can be found in the `NVCC Documentation <https://docs.nvidia.com/cuda/cuda-compiler-driver-nvcc/index.html#options-for-steering-gpu-code-generation>`_, 
along with details of the supported `virtual architectures <https://docs.nvidia.com/cuda/cuda-compiler-driver-nvcc/index.html#virtual-architecture-feature-list>`_ and `real architectures <https://docs.nvidia.com/cuda/cuda-compiler-driver-nvcc/index.html#gpu-feature-list>`_.

.. note:: SM 60 for Pascal GPUs is only available for CUDA 8.0 and above.

.. note:: SM 70 for Volta GPUs is only available for CUDA 9.0 and above.

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

   NVRM version: NVIDIA UNIX x86_64 Kernel Module  440.64.00  Wed Feb 26 16:26:08 UTC 2020
   GCC version:  gcc version 4.8.5 20150623 (Red Hat 4.8.5-39) (GCC)

Installation notes
------------------

These are primarily for system administrators.

Device driver
^^^^^^^^^^^^^

The NVIDIA device driver is installed and configured using the ``gpu-nvidia-driver`` systemd service (managed by Puppet).
This service runs ``/usr/local/scripts/gpu-nvidia-driver.sh`` at boot time to:

- Check the device driver version and uninstall it then reinstall the target version if required;
- Load the ``nvidia`` kernel module;
- Create several *device nodes* in ``/dev/``.

CUDA 10.2.89
^^^^^^^^^^^^

#. Installed with :download:`install.sh </sharc/software/install_scripts/libs/CUDA/install.sh>` with ``10.2.89_440.33.01`` as the sole argument. 
#. :download:`Modulefile </sharc/software/modulefiles/libs/CUDA/10.2.89/binary` was installed as ``/usr/local/modulefiles/libs/CUDA/10.2.89/binary``

CUDA 10.1.243
^^^^^^^^^^^^^

#. Installed with :download:`install.sh </sharc/software/install_scripts/libs/CUDA/install.sh>` with ``10.1.243_418.87.00`` as the sole argument. 
#. :download:`Modulefile </sharc/software/modulefiles/libs/CUDA/10.1.243/binary` was installed as ``/usr/local/modulefiles/libs/CUDA/10.1.243/binary``

CUDA 10.0.130
^^^^^^^^^^^^^

#. Installed with :download:`install.sh </sharc/software/install_scripts/libs/CUDA/install.sh>` with ``10.0.130_410.48`` as the sole argument. 
   This installs the toolkit and three NVIDIA-recommended patches.
#. :download:`Modulefile </sharc/software/modulefiles/libs/CUDA/10.0.130/binary>` was installed as ``/usr/local/modulefiles/libs/CUDA/10.0.130/binary``

CUDA 9.1.85
^^^^^^^^^^^

#. Installed with :download:`install.sh </sharc/software/install_scripts/libs/CUDA/install.sh>` with ``9.1.85_387.26`` as the sole argument. 
   This installs the toolkit and three NVIDIA-recommended patches.
#. :download:`Modulefile </sharc/software/modulefiles/libs/CUDA/9.1.85/binary>` was installed as ``/usr/local/modulefiles/libs/CUDA/9.1.85/binary``

CUDA 9.0.176
^^^^^^^^^^^^

#. Installed with :download:`install.sh </sharc/software/install_scripts/libs/CUDA/install.sh>` with ``9.0.176_384.81`` as the sole argument. 
   This installs the toolkit and four NVIDIA-recommended patches.
#. :download:`Modulefile </sharc/software/modulefiles/libs/CUDA/9.0.176/binary>` was installed as ``/usr/local/modulefiles/libs/CUDA/9.0.176/binary``

CUDA 8.0.44
^^^^^^^^^^^

#. Installed with :download:`install.sh </sharc/software/install_scripts/libs/CUDA/install.sh>` with ``8.0.44`` as the sole argument.
#. :download:`This modulefile </sharc/software/modulefiles/libs/CUDA/8.0.44/binary>` was installed as ``/usr/local/modulefiles/libs/CUDA/8.0.44/binary``

CUDA 7.5.18
^^^^^^^^^^^

#. Installed with :download:`install.sh </sharc/software/install_scripts/libs/CUDA/install.sh>` with ``7.5.18`` as the sole argument.
#. :download:`This modulefile </sharc/software/modulefiles/libs/CUDA/7.5.18/binary>` was installed as ``/usr/local/modulefiles/libs/CUDA/7.5.18/binary``
