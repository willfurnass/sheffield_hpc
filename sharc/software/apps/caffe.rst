.. _caffe_sharc:

Caffe
=====

.. sidebar:: Caffe

   :URL: http://caffe.berkeleyvision.org/

`Caffe <http://caffe.berkeleyvision.org/>`_ is a deep learning framework made with expression, speed, and modularity in mind. It is developed by the Berkeley Vision and Learning Center (`BVLC <http://bvlc.eecs.berkeley.edu/>`_) and by community contributors.

**Additional permissions are needed to use GPUs on Iceberg/ShARC. See** :ref:`GPUComputing_sharc` **for more information.**

Caffe Singularity Images
------------------------

Singularity images are self-contained virtual machines similar to Docker. For more information on Singularity and how to use the images, see :ref:`singularity_sharc`.

The following Singularity images are available on ShARC and can also be downloaded for use on your local machine:

* CPU Caffe rc5, Ubuntu 16.04, GCC 5.4.0
    * Image path: ``/usr/local/packages/singularity/caffe/rc5-CPU-Ubuntu16.04.img``
    * Def file: `Caffe CPU </sharc/software/apps/singularity/caffe_cpu.def>`
* GPU Caffe rc5, Ubuntu 16.04, CUDA 8, cuDNN 5.0, GCC 5.4.0
    * Image path: ``/usr/local/packages/singularity/caffe/rc5-GPU-Ubuntu16.04-CUDA8-cudNN5.0.img``
    * Def file: `Caffe GPU </sharc/software/apps/singularity/caffe_gpu.def>`

To get a bash terminal in to an image for example, use the command: ::

  singularity exec /usr/local/packages/singularity/caffe/rc5-CPU-Ubuntu16.04.img /bin/bash

The ``exec`` command can also be used to call any command/script inside the image e.g. ::

  singularity exec /usr/local/packages/singularity/caffe/rc5-CPU-Ubuntu16.04.img "caffe train your_solver.prototxt"

The paths ``/fastdata``, ``/data``, ``/home``, ``/scratch``, ``/shared`` are automatically mounted to your ShARC filestore directories. For GPU-enabled images the ``/nvlib`` and ``/nvbin`` is mounted to the correct Nvidia driver version for the node that you're using.

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

	module load module load apps/caffe/rc5/gcc-4.9.4-cuda-8.0-cudnn-5.1
	source activate caffe

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
