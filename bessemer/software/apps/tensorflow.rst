.. _tensorflow_bessemer:

Tensorflow
==========

.. sidebar:: Tensorflow

   :URL: https://www.tensorflow.org/

TensorFlow is an open source software library for numerical computation using data flow graphs. Nodes in the graph represent mathematical operations, while the graph edges represent the multidimensional data arrays (tensors) communicated between them. The flexible architecture allows you to deploy computation to one or more CPUs or GPUs in a desktop, server, or mobile device with a single API. TensorFlow was originally developed by researchers and engineers working on the Google Brain Team within Google's Machine Intelligence research organization for the purposes of conducting machine learning and deep neural networks research, but the system is general enough to be applicable in a wide variety of other domains as well.

About Tensorflow on Bessemer
----------------------------

**A GPU-enabled worker node must be requested in order to use the GPU version of this software. See** :ref:`GPUComputing_bessemer` **for more information.**

As Tensorflow and all its dependencies are written in Python, it can be installed locally in your home directory. The use of Anaconda (:ref:`python_conda_bessemer`) is recommended as it is able to create a virtual environment in your home directory, allowing the installation of new Python packages without admin permission.

This software and documentation is maintained by the `RSES group <https://rse.shef.ac.uk/>`_ and `GPUComputing@Sheffield <http://gpucomputing.shef.ac.uk/>`_. For feature requests or if you encounter any problems, please raise an issue on the `GPU Computing repository <https://github.com/RSE-Sheffield/GPUComputing/issues>`_.



Installation in Home Directory - CPU Version
--------------------------------------------

In order to to install to your home directory, Conda is used to create a virtual python environment for installing your local version of Tensorflow.

First request an interactive session, e.g. with :ref:`slurm_interactive`.

Then Tensorflow can be installed by the following ::

  #Load the conda module
  module load Anaconda3/5.3.0 

  #Create an conda virtual environment called 'tensorflow'
  conda create -n tensorflow python=3.6

  #Activate the 'tensorflow' environment
  source activate tensorflow

  pip install tensorflow


**Every Session Afterwards and in Your Job Scripts**

Every time you use a new session or within your job scripts, the modules must be loaded and conda must be activated again. Use the following command to activate the Conda environment with Tensorflow installed: ::

  module load Anaconda3/5.3.0 
  source activate tensorflow


Installation in Home Directory - GPU Version
--------------------------------------------

The GPU version of Tensorflow comes in a different PIP package and is also dependent on CUDA and cuDNN libraries making the installation procedure slightly different.

First request an interactive session, e.g. see :ref:`GPUInteractive_bessemer`.

Then GPU version of Tensorflow can be installed by the following ::

  #Load the conda module
  module load Anaconda3/5.3.0 

  #Load the CUDA and cuDNN module
  module load cuDNN/7.4.2.24-gcccuda-2019a

  #Create an conda virtual environment called 'tensorflow-gpu'
  conda create -n tensorflow-gpu python=3.6

  #Activate the 'tensorflow-gpu' environment
  source activate tensorflow-gpu

  #Install GPU version of Tensorflow 
  pip install tensorflow-gpu

If you wish to use an older version of tensorflow-gpu, you can do so using :code:`pip install tensorflow-gpu==<version_number>`

**Every Session Afterwards and in Your Job Scripts**

Every time you use a new session or within your job scripts, the modules must be loaded and conda must be activated again. Use the following command to activate the Conda environment with Tensorflow installed: ::

  module load Anaconda3/5.3.0 
  module load cuDNN/7.4.2.24-gcccuda-2019a
  source activate tensorflow-gpu


Testing your Tensorflow installation
------------------------------------

You can test that Tensorflow is running on the GPU with the following python code ::

  import tensorflow as tf
  # Creates a graph
  #If using CPU, replace /device:GPU:0 with /cpu:0
  with tf.device('/device:GPU:0'):
    a = tf.constant([1.0, 2.0, 3.0, 4.0, 5.0, 6.0], shape=[2, 3], name='a')
    b = tf.constant([1.0, 2.0, 3.0, 4.0, 5.0, 6.0], shape=[3, 2], name='b')
    c = tf.matmul(a, b)
  # Creates a session with log_device_placement set to True.
  sess = tf.Session(config=tf.ConfigProto(log_device_placement=True))
  # Runs the op.
  print(sess.run(c))

Which gives the following results ::

	[[ 22.  28.]
	 [ 49.  64.]]

CUDA and CUDNN Import Errors
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Tensorflow releases depend on specific versions of both CUDA and CUDNN. If the wrong CUDNN module is loaded, you may receive an :code:`ImportError` runtime errors such as: 

.. code-block :: python

   ImportError: libcublas.so.10.0: cannot open shared object file: No such file or directory


This indicates that Tensorflow was expecting to find CUDA 10.0 (and an appropriate version of CUDNN) but was unable to do so.

The following table shows the which module to load for the various versions of Tensorflow, based on the `tested build configurations <https://www.tensorflow.org/install/source#linux>`_. 



+------------+------+--------+--------------------------------------------+
| Tensorflow | CUDA | CUDNN  | Module                                     | 
+============+======+========+============================================+
| 2.1.0      | 10.1 | >= 7.4 | `libs/cudnn/7.6.5.32/binary-cuda-10.1.243` |
+------------+------+--------+--------------------------------------------+
| 2.0.0      | 10.0 | >= 7.4 | `libs/cudnn/7.5.0.56/binary-cuda-10.0.130` |
+------------+------+--------+--------------------------------------------+
| 1.14.0     | 10.0 | >= 7.4 | `libs/cudnn/7.5.0.56/binary-cuda-10.0.130` |
+------------+------+--------+--------------------------------------------+
| 1.13.1     | 10.0 | >= 7.4 | `libs/cudnn/7.5.0.56/binary-cuda-10.0.130` |
+------------+------+--------+--------------------------------------------+
| >= 1.5.0   |  9.0 | 7      | `libs/cudnn/7.3.1.20/binary-cuda-9.0.176`  |
+------------+------+--------+--------------------------------------------+
| >= 1.3.0   |  8.0 | 6      | `libs/cudnn/6.0/binary-cuda-8.0.44`        |
+------------+------+--------+--------------------------------------------+
| >= 1.0.0   |  8.0 | 5.1    | `libs/cudnn/5.1/binary-cuda-8.0.44`        |
+------------+------+--------+--------------------------------------------+



