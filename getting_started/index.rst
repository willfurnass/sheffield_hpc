.. _getting-started:


Getting Started
===============

If you have not used a High Performance Computing (HPC) cluster, Linux or
even a command line before this is the place to start. This guide will get you 
set up using iceberg in the easiest way that fits your requirements.

Getting an Account
##################

Before you can start using iceberg you need to register for an account.
Accounts are availible for staff by emailing `helpdesk@sheffield.ac.uk <helpdesk@sheffield.ac.uk>`_.

The following categories of students can also have an account on iceberg with 
the permission of their supervisors:

* Research Postgraduates
* Taught Postgraduates - project work
* Undergraduates 3rd & 4th year  - project work

Student's supervisors can request an account for the student by emailing
`helpdesk@sheffield.ac.uk <helpdesk@sheffield.ac.uk>`_.

.. note::

    Once you have obtained your iceberg username, you need to initialize your 
    iceberg password to be the same as your normal password by using the CICS
    `synchronize passwords system <https://www.shef.ac.uk/cics/password>`_.

What do you want to Achieve?
############################

The iceberg computing cluster is a very powerful and flexible service. 
It can be used for running many different types of computationally expensive 
tasks. The best way for you to use iceberg, and the tools you will need to learn
depend on what you want to do.




Q) Does your task take less than 8 hours to run?

   If you have a 'short' running task, which takes less than 8 hours, either 
   a graphical application such as MATLAB or Ansys or a script, you can 
   use an interactive session. See :ref:`short-running`.

Q) Does your task take longer than 8 hours, or do you have lots of short tasks?

   If you need more than 8 hours of compute time, or you need to run lots of 
   small tasks you will need to learn the Grid Engine submission system.
   See :ref:`sge-intro`.

Q) Does your graphical or interactive taks require multiple threads or lots of memory?

   You will need to specify Grid Engine parameters to your interactive job 
   command. See :ref:`sge-interactive`.

Q) Do you need to run a parallel application on multiple processors or nodes?

   If you have a task that uses distributed memory parallelism, like ePI, you
   will need to read the section on ``mpie``.

Q) Do you want to visualise your results without downloading all the data from iceberg?

   You can make use of iceberg's remote visualisation nodes, or run 
   paraview or the Jupyter (IPython) notebook in server mode. See ``remote-vis``.


.. _short-running:

Short Running Applications
##########################

.. _connecting:

Terminal Access
```````````````

Accessing iceberg through a terminal is the most flexible way of using iceberg,
as it is the native way of interfacing with the linux cluster.


Open a Terminal
---------------

**Windows**


Download and install `mobaxterm <https://mobaxterm.mobatek.net/download-home-edition.html>`_.

Running MobaXterm should display the following screen:

.. image:: /images/mobaxterm-welcome.png
   :width: 50%
   :align: center

Once you have this screen you can jump to :ref:`ssh`.


**Mac OS/X and Linux**

Linux and Mac OS/X both have a terminal emulator program pre-installed.

Open a terminal and then go to :ref:`ssh`.

.. _ssh:

Connect to iceberg
------------------

Once you have a terminal open run the following command::

    ssh -X <username>@iceberg.shef.ac.uk

where you replace `<username>` with your CICS username.

This should give you a prompt resembling the one below::

    [te1st@iceberg-login2 ~]$ 

at this prompt type::

    qsh

like this::

    [te1st@iceberg-login2 ~]$ qsh
    Your job 135355 ("INTERACTIVE") has been submitted
    waiting for interactive job to be scheduled ....
    Your interactive job 135355 has been successfully scheduled.

which will pop up another terminal window, which supports graphical applications.

.. note::

    Iceberg is a compute cluster. When you login to the cluster you reach one 
    of two login nodes. You **should not** run applications on the login nodes.
    Running qsh gives you an interactive terminal on one of the many worker nodes
    in the cluster.

    If you only need terminal based (CLI) applications you can run the qrsh command.
    Which will give you a shell on a worker node, but without graphical application
    (X server) support.


Now go to :ref:`running-applications`.

.. raw:: html

   <p>
    This video shows the connection process using mobaxterm, and then connection 
    and running matlab from a <cite>qsh</cite> terminal.
   </p>

   <video style="margin-left: auto; margin-right:auto; display: block;" width=70% controls>
       <source src="http://rcg.group.shef.ac.uk/tutorial_videos/mobaxterm-login-matlab-demo.webm" type="video/webm" />
       <source src="http://rcg.group.shef.ac.uk/tutorial_videos/mobaxterm-login-matlab-demo.mp4" type="video/mp4" />
   </video>


#Option 2: Web Interface
#```````````````````````
#
#Not me
#
#


.. _running-applications:

Loading Applications or Libraries
---------------------------------

Once you have connected to iceberg you have to run `qsh` or `qrsh` to get a 
terminal on a woker node (see :ref:`ssh`).
Once you have a terminal on a worker node you can load software using the modules
system.
For information on modules and the applications availible see :ref:`software`.



.. toctree::
    :hidden:

    self
    sge-intro
