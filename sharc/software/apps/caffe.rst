.. _caffe_sharc:

Caffe
=====

.. sidebar:: Caffe

   :URL: http://caffe.berkeleyvision.org/

`Caffe <http://caffe.berkeleyvision.org/>`_ is a deep learning framework made with expression, speed, and modularity in mind. It is developed by the Berkeley Vision and Learning Center (`BVLC <http://bvlc.eecs.berkeley.edu/>`_) and by community contributors.

About Caffe on ShARC
--------------------

Caffe is available on ShARC as both Singularity images and as a module.

This software and documentation is maintained by the `RSES group <http://rse.shef.ac.uk/>`_ and `GPUComputing@Sheffield <http://gpucomputing.shef.ac.uk/>`_. For feature requests or if you encounter any problems, please raise an issue on the `GPU Computing repository <https://github.com/RSE-Sheffield/GPUComputing/issues>`_.


Caffe Singularity Images
------------------------

Singularity images are self-contained virtual machines similar to Docker. For more information on Singularity and how to use the images, see :ref:`singularity_sharc`.

A symlinked file is provided that always point to the latest image: ::

  #CPU Caffe
  /usr/local/packages/singularity/images/caffe/cpu.img

  #GPU Caffe
  /usr/local/packages/singularity/images/caffe/gpu.img

To get a bash terminal in to an image for example, use the command: ::

  singularity exec /usr/local/packages/singularity/images/caffe/gpu.img /bin/bash

The ``exec`` command can also be used to call any command/script inside the image e.g. ::

  singularity exec /usr/local/packages/singularity/images/caffe/gpu.img caffe train -solver=your_solver.prototxt

The paths ``/fastdata``, ``/data``, ``/home``, ``/scratch``, ``/shared`` are automatically mounted to your ShARC filestore directories. For GPU-enabled images the ``/nvlib`` and ``/nvbin`` is mounted to the correct Nvidia driver version for the node that you're using.

Image Index
^^^^^^^^^^^

Paths to the actual images and definition files are provided below for downloading and building of custom images.

* Shortcut to Latest Image
    * CPU
        * ``/usr/local/packages/singularity/images/caffe/cpu.img``
    * GPU
        * ``/usr/local/packages/singularity/images/caffe/gpu.img``
* Images list
    * CPU
        * Latest: rc5-CPU-Ubuntu16.04 (Python 2.7)
            * Path: ``/usr/local/packages/singularity/images/caffe/rc5-CPU-Ubuntu16.04.img``
            * Def file: `/sharc/software/apps/singularity/caffe_cpu.def </sharc/software/apps/singularity/caffe_cpu.def>`
    * GPU
        * Latest: rc5-GPU-Ubuntu16.04-CUDA8-cudNN5.0 (Python 2.7)
            * Path: ``/usr/local/packages/singularity/images/caffe/rc5-GPU-Ubuntu16.04-CUDA8-cudNN5.0.img``
            * Def file: `/sharc/software/apps/singularity/caffe_gpu.def </sharc/software/apps/singularity/caffe_gpu.def>`

Using the Caffe Module
----------------------

First request a GPU interactive session (see :ref:`GPUInteractive_sharc`).

The Caffe module can be loaded with the following command:   ::

  module load apps/caffe/rc5/gcc-4.9.4-cuda-8.0-cudnn-5.1

Installing Additional Python Modules
------------------------------------

The Caffe module is pre-installed with Anaconda version 3.4.2. You can install additional python packages by creating a virtual python environment in your home directory using conda. ::

  #Creates a conda environment named caffe
	conda create -n caffe python=3.5
  #Activates the caffe python environment
	source activate caffe

You will also need to install ``numpy`` which can be obtained from the conda repository. ::

	conda install -c anaconda numpy=1.11.2


Every Session Afterwards and in Your Job Scripts
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
If you created a virtual python environment, you must activate it at every new session and within your job scripts: ::

	module load apps/caffe/rc5/gcc-4.9.4-cuda-8.0-cudnn-5.1
	source activate caffe

Caffe Training
--------------

`GPUComputing@sheffield <http://gpucomputing.shef.ac.uk>`_ provides training materials on the `use of Caffe on the DGX-1 and ShARC cluster <http://gpucomputing.shef.ac.uk/education/intro_dl_sharc_dgx1/>`_.

Installation Notes
------------------

For the module: ::

  module load apps/caffe/rc5/gcc-4.9.4-cuda-8.0-cudnn-5.1

The following modules are automatically loaded:
  * GCC 4.9.4
  * CUDA 8
  * cuDNN 5.1

And comes with the following libraries:
  * Anaconda 4.2.0 (Python 3)
  * boost
  * protobuf
  * hdf5
  * snappy
  * glog
  * gflags
  * openblas
  * leveldb
  * lmdb
  * yasm
  * libx264
  * libx265
  * libfdk_acc
  * libopus
  * libogg
  * libvorbis
  * freetype
  * ffmpeg
  * libjpeg
  * libpng
  * libtiff
  * opencv 3.2.0
