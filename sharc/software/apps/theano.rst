.. _theano_sharc:

Theano
======

.. sidebar:: Theano

   :URL: http://deeplearning.net/software/theano/index.html

Theano is a Python library that allows you to define, optimize, and evaluate mathematical expressions involving multi-dimensional arrays efficiently. Theano is most commonly used to perform deep learning and has excellent GPU support and integration through PyCUDA. The following steps can be used to setup and configure Theano on your own profile.

About Theano on ShARC
---------------------

As Theano and all its dependencies are written in Python, it can be installed locally in your home directory. The use of Anaconda (:ref:`sharc-python-conda`) is recommended as it is able to create a virtual environment in your home directory, allowing the installation of new Python packages without admin permission.

**Additional permissions are needed to use GPUs on Iceberg/ShARC. See** :ref:`GPUIntro` **for more information.**


Installation
------------

First request an interactive session with a sufficient amount of memory, e.g. ::

		qrshx -l gpu=1 -l gpu_arch=nvidia-k40m -l mem=13G


Load the relevant modules with the following command: ::

	module load apps/python/anaconda3-4.2.0
	module load libs/cudnn/5.1/binary-cuda-8.0.44

Create a conda environment to load relevant modules on your local user account ::

		conda create -n theano python=3.5
		source activate theano

Install the other Python module dependencies which are required using pip (alternatively these could be installed with conda if you prefer) ::

		pip install theano nose nose-parameterized pycuda



For optimal Theano performance, enable the CUDA memory manager CNMeM. To do this, create the .theanorc file in your HOME directory and set the fraction of GPU memory reserved by Theano. The exact amount of energy may have to be hand-picked: if Theano asks for more memory that is currently available on the GPU, an error will be thrown during import of theano module. Create or edit the ``.theanorc`` file with nano ::

		nano ~/.theanorc

Add the following lines and, if necessary, change the 0.8 number to whatever works for you ::

		[lib]
		cnmem=0.8

Run python and verify that Theano is working correctly ::

		python -c "import theano;theano.test()"

Running Theano in Another Session or in the Job Script
------------------------------------------------------

The previous instuctions installs Theano and its dependencies inside your home directory but every time you use a new session or within your job scripts, the modules must be loaded and conda must be activated again. Use the following command to activate the Conda environment with Theano installed: ::

	module load apps/python/anaconda3-4.2.0
	module load libs/cudnn/5.1/binary-cuda-8.0.44
	source activate theano
