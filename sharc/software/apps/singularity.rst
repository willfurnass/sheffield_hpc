.. _singularity_sharc:

Singularity
===========

.. sidebar:: Singularity

   :Version: 2.3.2
   :URL: http://singularity.lbl.gov/

Designed around the notion of extreme mobility of compute and reproducible science,
Singularity enables users to have full control of their operating system environment.
This means that a non-privileged user can "swap out" the operating system on the host for one they control.
So if the host system is running Centos Linux but your application runs in Ubuntu Linux,
you can create an Ubuntu image,
install your applications into that image,
copy the image to another host,
and run your application on that host in it’s native Ubuntu environment.

Singularity also allows you to leverage the resources of whatever host you are on.
This includes high-speed cluster interconnects,
resource managers,
file systems,
GPUs and/or
accelerators, etc.

About Singularity Containers (Images)
-------------------------------------

Similar to Docker,
a Singularity container (image) is a self-contained software stack.
As Singularity does not require a root-level daemon to run its images
it is compatible for use with ShARC's scheduler inside your job scripts.
The running images also uses the credentials of the person calling it.

In practice, this means that an image created on your local machine
with all your research software installed for local development
will also run on the ShARC cluster.

Pre-built images have been provided on the cluster and
can also be download for use on your local machine
(see :ref:`use_image_singularity_sharc`).
Creating and modifying images however,
requires root permission and so
must be done on your machine (see :ref:`create_image_singularity_sharc`).

Singularity images are currently provided for:

* :ref:`caffe_sharc`
* :ref:`theano_sharc`
* :ref:`torch_sharc`
* :ref:`tensorflow_sharc`

.. _use_image_singularity_sharc:

Interactive Usage of Singularity Images
---------------------------------------

**To use Singularity interactively, an interactive session must first be requested using** :ref:`qrshx` **for example.**

To get an interactive shell in to the image, use the following command: ::

  singularity shell path/to/imgfile.img

Or if you prefer bash: ::

  singularity exec path/to/imgfile.img /bin/bash

Note that the ``exec`` command can also be used to execute other applications/scripts inside the image or
from the mounted directories (See :ref:`auto_mounting_filestore_singularity_sharc`): ::

    singularity exec path/to/imgfile.img my_script.sh

.. note::

    You may get a warning similar to:

    .. code-block:: none

        groups: cannot find name for group ID ...

    :ref:`This can be ignored <unnamed_groups>` and will not have an affect on running the image.

.. _use_image_batch_singularity_sharc:

Submitting Batch Jobs That Uses Singularity Images
--------------------------------------------------

When submitting a job that uses a Singularity image,
it is not possible to use the interactive shell
(e.g. ``singularity shell`` or ``singularity exec path/to/imgfile.img /bin/bash``).
You must use the ``exec`` command to call the desired application or script directly.

For example, if we would like to use a command ``ls /`` to get the content of the root folder in the image,
two approaches are shown in the following job script ``my_singularity_job.sh``:

.. code-block:: bash

  #!/bin/bash
  #$ -l rmem=8G
  # We requested 8GB of memory in the line above, change this according to your
  # needs e.g. add -l gpu=1 to reqest a single GPU

  #Calling ls directly using the exec command
  singularity exec path/to/imgfile.img ls /

  #Have Singularity call a custom script from your home or other mounted directories
  #Don't forget to make the script executable before running by using chmod
  chmod +x ~/myscript.sh
  singularity exec path/to/imgfile.img ~/myscript.sh

Where the content of ``~/myscript.sh`` is shown below:

.. code-block:: bash

  #!/bin/bash

  ls /

The job can then be submitted as normal with ``qsub``: ::

  qsub my_singularity_job.sh

.. _auto_mounting_filestore_singularity_sharc:

Automatic Mounting of ShARC Filestore Inside Images
----------------------------------------------------

When running Singularity images on ShARC,
the paths ``/fastdata``, ``/data``, ``/home``, ``/scratch``, ``/shared`` are
automatically mounted to your ShARC directories.

Images that uses the GPU requires driver files that matches the host system.
In ShARC these files are located outside of the image and
automatically mounted to paths ``/nvbin`` and ``/nvlib`` within the image.

Image Index on Github
---------------------

All Singularity container definitions available on ShARC can be found at `https://github.com/rses-singularity <https://github.com/rses-singularity>`_. The definition files can be used as a template for building your own images.


Installing Singularity on Your Local Machine
--------------------------------------------

You will need Singularity installed on your machine in order to
locally run, create and modify images.
The following is the installation command for debian/ubuntu based systems:

.. code-block:: bash

  VERSION=2.3.1
  wget https://github.com/singularityware/singularity/releases/download/$VERSION/singularity-$VERSION.tar.gz
  tar xvf singularity-$VERSION.tar.gz
  cd singularity-$VERSION
  ./configure --prefix=/usr/local
  make
  sudo make install


Manually mounting paths
-----------------------

When using ShARC's pre-built images on your local machine,
it may be useful to mount the existing directories in the image to your own path.
This can be done with the flag ``-B local/path:image/path`` with
the path outside of the image left of the colon and
the path in the image on the right side, e.g. ::

  singularity shell -B local/datapath:/data,local/fastdatapath:/fastdata path/to/imgfile.img

The command mounts the path ``local/datapath`` on your local machine to
the ``/data`` path in the image.
Multiple mount points can be joined with ``,``
as shown above where we additionally specify that ``local/fastdata`` mounts to ``/fastdata``.
The ``/home`` folder is automatically mounted by default.

**Note: In order to mount a path, the directory must already exist within the image.**

.. _create_image_singularity_sharc:

Creating Your Own Singularity Images
------------------------------------

**Root access is required for modifying Singularity images so if you need to edit an image it must be done on your local machine.  However you can create disk images and import docker images using normal user privileges on recent versions of Singularity**

Firstly an empty image must be created. The following command creates an image named ``myimage.img`` of the size 1024 MB: ::

  singularity create -s 1024 myimage.img

Singularity uses a definition file for bootstrapping an image. An example definition ``ShARC-Ubuntu-Base.def`` is shown below ::

  Bootstrap: docker
  From: ubuntu:latest

  %setup
  	#Runs on host. The path to the image is $SINGULARITY_ROOTFS

  %post
  	#Post setup, runs inside the image

    #Default mount paths
  	mkdir /scratch /data /shared /fastdata

    #Nvidia driver mount paths, only needed if using GPU
  	mkdir /nvlib /nvbin

    #Add nvidia driver paths to the environment variables
  	echo "\n #Nvidia driver paths \n" >> /environment
  	echo 'export PATH="/nvbin:$PATH"' >> /environment
  	echo 'export LD_LIBRARY_PATH="/nvlib:$LD_LIBRARY_PATH"' >> /environment

  %runscript
    #Runs inside the image every time it starts up

  %test
    #Test script to verify that the image is built and running correctly

The definition file takes a base image from `docker hub <https://hub.docker.com/>`_,
in this case the latest version of Ubuntu ``ubuntu:latest``.
Other images on the hub can also be used as the base for the Singularity image,
e.g. ``From: nvidia/cuda:8.0-cudnn5-devel-ubuntu16.04`` uses Nvidia's docker image with Ubuntu 16.04 that already has CUDA 8 installed.

After creating a definition file, use the ``bootstrap`` command to build the image you've just created: ::

  sudo singularity bootstrap myimage.img ShARC-Ubuntu-Base.def

You can also modify the contents of an image after it's been created using the ``-w`` flag: ::

  sudo singularity shell -w myimage.img

The command above gives you a shell in to the image with root access that can then be used to modify its contents.

If you just want to import a docker image directly without making any modifications to it you can run the following command without requiring sudo privileges: ::

  singularity import myimage.img docker://ubuntu:latest

Using Nvidia GPU with Singularity Images on Your Local Machine
--------------------------------------------------------------

**Support is only available for machines with Nvdia GPUs and will not work for other GPU manufacturers (e.g. AMD).**

In order to use Nvidia GPUs within a singularity image,
a copy of the driver files must be present in the image and
must match the version of the host machine.
`Previously <https://hpc.nih.gov/apps/singularity.html>`_,
this is done by embedding the driver within the image itself
which creates a non-portable image.

On the ShARC cluster,
these driver files are stored outside of the image and
automatically mounted to the folders ``/nvbin`` and ``/nvlib`` at run-time.
To use the images locally on your machine you simply need to
provide the correct driver files for the machine you're using.

Use the following command to find your current driver version: ::

  nvidia-smi

Where you will get something similar to the following:

.. code-block:: none

  Tue Mar 28 16:43:08 2017
  +-----------------------------------------------------------------------------+
  | NVIDIA-SMI 367.57                 Driver Version: 367.57                    |
  |-------------------------------+----------------------+----------------------+
  | GPU  Name        Persistence-M| Bus-Id        Disp.A | Volatile Uncorr. ECC |
  | Fan  Temp  Perf  Pwr:Usage/Cap|         Memory-Usage | GPU-Util  Compute M. |
  |===============================+======================+======================|
  |   0  GeForce GTX TITAN   Off  | 0000:01:00.0      On |                  N/A |
  | 30%   35C    P8    18W / 250W |    635MiB /  6078MiB |      1%      Default |
  +-------------------------------+----------------------+----------------------+

It can be seen that the driver version on our current machine is ``367.57``.
Go to the `Nvidia website <http://nvidia.com>`_ and
search for the correct Linux driver for your graphics card.
Download the :download:`extract_nvdriver_and_moveto.sh </sharc/software/apps/singularity/extract_nvdriver_and_moveto.sh>` to
the same directory and run it like so: ::

  chmod +x extract_nvdriver_and_moveto.sh
  extract_driver_and_moveto.sh 367.57 ~/mynvdriver

If you're using the Singularity definition file
as shown above (see :ref:`create_image_singularity_sharc`),
the ``/nvbin`` and ``/nvlib`` directories will have been created.
They simply need to be correctly mounted when
running the image using the command where
our extracted driver files are located at ``~/mynvdriver``: ::

  singularity shell -B ~/mynvdriver:/nvlib,~/mynvdriver:/nvbin myimage.img

How Singularity is installed and 'versioned' on the cluster
-----------------------------------------------------------

Singularity, unlike much of the other key software packages on ShARC,
is not activated using module files.
This is because module files are primarily for the purpose of
being able to install multiple version of the same software
and for security reasons only the most recent version of Singularity is installed.
The security risks associated with providing outdated builds of Singularity
are considered to outweigh the risk of upgrading to backwards incompatible versions.

Singularity has been installed on all worker nodes
using the latest RPM package
from the `EPEL <https://fedoraproject.org/wiki/EPEL>`_ repository.
