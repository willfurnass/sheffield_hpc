.. _caffe_sharc:

Caffe
=====

.. sidebar:: Caffe

   :URL: http://caffe.berkeleyvision.org/

`Caffe <http://caffe.berkeleyvision.org/>`_ is a deep learning framework made with expression, speed, and modularity in mind. It is developed by the Berkeley Vision and Learning Center (`BVLC <http://bvlc.eecs.berkeley.edu/>`_) and by community contributors.

**Additional permissions are needed to use GPUs on Iceberg/ShARC. See** :ref:`GPUComputing_sharc` **for more information.**

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
