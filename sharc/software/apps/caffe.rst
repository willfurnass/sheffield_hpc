.. _caffe_sharc:

Caffe
=====

.. sidebar:: Caffe

   :URL: http://caffe.berkeleyvision.org/

`Caffe <http://caffe.berkeleyvision.org/>`_ is a deep learning framework made with expression, speed, and modularity in mind. It is developed by the Berkeley Vision and Learning Center (`BVLC <http://bvlc.eecs.berkeley.edu/>`_) and by community contributors.

About Caffe on ShARC
--------------------



Using the Caffe Module
----------------------

There are two version of Caffe on ShARC (see documentation), the standard version that is newer (has more layer types) and Nvidia’s version that is more optimised for multi-GPU use. We’ll be using the standard version as it has support for recurrent layers.

To use the python bindings for Caffe we need to prepare our python environment. First we load up the Caffe module which also loads Conda 3.4, CUDA 8, cuDNN 5.1 and GCC 4.9.4 ::

  module load libs/caffe/nvidia-0.15.13/gcc-4.9.4-cuda-8.0-cudnn-5.1-python-2.7-conda-3.4-TESTING
  module load libs/caffe/rc3/gcc-4.9.4-cuda-8.0-cudnn-5.1-conda-3.4-TESTING
