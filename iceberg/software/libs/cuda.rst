.. _`cuda_iceberg`:

CUDA
====
CUDA, which stands for Compute Unified Device Architecture, is a parallel computing platform and application programming interface (API) model created by NVIDIA.
It allows software developers to use a CUDA-enabled graphics processing unit (GPU) for general purpose processing - an approach known as *General Purpose GPU* (GPGPU).

Usage
-----
There are several versions of the CUDA library available. As with many libraries installed on the system, CUDA libraries are made available via ``module`` commands which are only available once you have started a ``qrsh`` or ``qsh`` session.

The latest version CUDA is loaded with the command ::

        module load libs/cuda

Alternatively, you can load a specific version with one of the following: ::

        module load libs/cuda/8.0.44
        module load libs/cuda/7.5.18
        module load libs/cuda/6.5.14
        module load libs/cuda/4.0.17
        module load libs/cuda/3.2.16

To check which version of CUDA you are using ::

        $ nvcc --version
        nvcc: NVIDIA (R) Cuda compiler driver
        Copyright (c) 2005-2016 NVIDIA Corporation
        Built on Sun_Sep__4_22:14:01_CDT_2016
        Cuda compilation tools, release 8.0, V8.0.44

**Important** To compile CUDA programs you also need a compatible version of the :ref:`<gcc_iceberg>`.  As of version 8.0.44, CUDA is compatible with GCC versions:

* greater than or equal to 4.7.0 (to allow for the use of c++11 features) and
* less than 5.0.0

It is therefore recommended that you load the most recent 4.x version of GCC when building CUDA programs on Iceberg: ::

        module load compilers/gcc/4.9.2

Compiling the sample programs
-----------------------------
You do not need to be using a GPU-enabled node to compile the sample programs but you do need a GPU to run them.

In a `qrsh` session ::

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

A basic test is to run one of the resulting binaries, ``deviceQuery``.

Documentation
-------------
* `CUDA Toolkit Documentation <http://docs.nvidia.com/cuda/index.html#axzz3uLoSltnh>`_
* `The power of C++11 in CUDA 7 <http://devblogs.nvidia.com/parallelforall/power-cpp11-cuda-7/>`_

CUDA Training
-------------

`GPUComputing@sheffield <http://gpucomputing.shef.ac.uk>`_ provides a self-paced `introduction to CUDA <http://gpucomputing.shef.ac.uk/education/cuda/>`_ training course.

Determining the NVIDIA Driver version
-------------------------------------
Run the command: ::

        cat /proc/driver/nvidia/version

Example output is ::

        NVRM version: NVIDIA UNIX x86_64 Kernel Module  367.44  Wed Aug 17 22:24:07 PDT 2016
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

The NVIDIA device driver is currently version 367.44.  The driver installer provides OpenGL libraries.

CUDA 8.0.44
^^^^^^^^^^^

#. The CUDA toolkit binaries and samples were installed using a binary ``.run`` file: ::

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
