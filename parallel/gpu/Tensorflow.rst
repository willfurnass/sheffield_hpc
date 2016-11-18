.. _Tensorflow:

Deep Learning with Tensorflow on GPU nodes
------------------------------------------

TensorFlow is an open source software library for numerical computation using data flow graphs. Nodes in the graph represent mathematical operations, while the graph edges represent the multidimensional data arrays (tensors) communicated between them. The flexible architecture allows you to deploy computation to one or more CPUs or GPUs in a desktop, server, or mobile device with a single API. TensorFlow was originally developed by researchers and engineers working on the Google Brain Team within Google's Machine Intelligence research organization for the purposes of conducting machine learning and deep neural networks research, but the system is general enough to be applicable in a wide variety of other domains as well.

The following is an instruction on how to setup Tensorflow on your local user account


Request an interactive session with a sufficient amount of memory ::

		qsh -l gpu=1 --l gpu_arch=nvidia-k40m -l mem=13G

Load the relevant modules (our example uses CUDA 8.0 with cuDNN 5.0 but :ref:`other versions are available <iceberg_cudnn>`) ::

		module load apps/python/anaconda3-2.5.0
		module load libs/binlibs/cudnn/cuda8.0/cudnn5.0
		module load compilers/gcc/5.4
		module load apps/java/1.8u71


Create a conda environment to load relevant modules on your local user account ::

		conda create -n tensorflow python=3.5 anaconda3-2.5.0 
		source activate tensorflow
		
What next??...