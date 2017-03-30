.. _getting-started:

Getting Started
===============

If you have not used a High Performance Computing (HPC) cluster, the Linux operating system (an alternative to Windows and macOS) or
even a command line before this is the place to start. This guide will get you 
set up using the University's clusters in the easiest way that fits your
requirements.

Getting an Account
##################

Before you can start the clusters you need to register for an account.
Accounts are availible for staff by emailing `helpdesk@sheffield.ac.uk <helpdesk@sheffield.ac.uk>`_.

The following categories of students can also have an account with 
the permission of their supervisors:

* Research Postgraduates
* Taught Postgraduates - project work
* Undergraduates 3rd & 4th year  - project work

Student's supervisors can request an account for the student by emailing
`helpdesk@sheffield.ac.uk <helpdesk@sheffield.ac.uk>`_.

.. note::

    Once you have obtained your username, you need to initialize your 
    password to be the same as your normal password by using the CICS
    `synchronize passwords system <https://www.shef.ac.uk/cics/password>`_.

.. _connecting:

Connecting to a cluster
#######################

The most versatile way to **run commands and submit jobs** on one of the clusters is to use a mechanism called `SSH <https://en.wikipedia.org/wiki/Secure_Shell>`__, 
which is a common way of remotely logging in to computers running the Linux operating system.

To connect to another machine using SSH you need to have a SSH *client* program installed on your machine.  
macOS and Linux come with a command-line (text-only) SSH client pre-installed.  
On Windows there are various graphical SSH clients you can used, including *MobaXTerm*.
See below for more information on how to use these clients to connect to our clusters.

To **transfer files to/from the clusters** you can either:

* Use a program that supports one or both of the SCP and SFTP *protocols* to copy/move files to/from your own machine
* Use a `Research Storage fileshare <http://www.sheffield.ac.uk/cics/research-storage/>`_ as common storage directly 
  accessible from your own machine and from the clusters.

A more comprehensive guide to accessing the clusters and transfering files is located at 
`http://www.sheffield.ac.uk/cics/research/hpc/using/access <http://www.sheffield.ac.uk/cics/research/hpc/using/access>`_

.. note::

    You can connect to **Iceberg** using SSH/SCP/SFTP from anywhere (on campus and off campus with/without a `VPN connection <https://www.sheffield.ac.uk/cics/vpn>`_)
    but if you are using **ShARC** and you are **off-campus** then you need to set up a VPN connection first.
    Instructions for setting up a VPN connection are `available for Windows, macOS and Linux <https://www.sheffield.ac.uk/cics/vpn>`_.

Connecting from Windows
```````````````````````

Download and install the *Installer edition* of `mobaxterm <https://mobaxterm.mobatek.net/download-home-edition.html>`_.

After starting MobaXterm you should see something like this:

.. image:: /images/mobaxterm-welcome.png
   :width: 50%
   :align: center

Click *Start local terminal* and if you see something like the the following then please continue to the :ref:`ssh` section.

.. image:: /images/mobaxterm-terminal.png
   :width: 50%
   :align: center


Mac OS/X and Linux
``````````````````

Linux and macOS (OS X) both have a terminal emulator program pre-installed.  
If you are using macOS and want to be able to run graphical applications on the clusters then 
you need to install the latest version of the `XQuartz <https://www.xquartz.org/>`_ *X Windows server*.

Open a terminal and then go to :ref:`ssh`.

.. _ssh:

Connect to iceberg
``````````````````

Once you have a terminal open run the following command: ::

    ssh -X <username>@iceberg.shef.ac.uk

where you replace `<username>` with your CICS username.

.. note::

    **macOS users**: if this fails then:
    
    * Check that your `XQuartz <https://www.xquartz.org/>`_ is up to date then try again *or*
    * Try again with ``-Y`` instead of ``-X``

This should give you a prompt resembling the one below: ::

    [te1st@iceberg-login1 ~]$ 

at this prompt type: ::

    qsh

like this: ::

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


.. raw:: html

   <p>
    This video shows the connection process using mobaxterm, and then connection 
    and running matlab from a <cite>qsh</cite> terminal.
   </p>

   <video style="margin-left: auto; margin-right:auto; display: block;" width=70% controls>
       <source src="http://rcg.group.shef.ac.uk/tutorial_videos/mobaxterm-login-matlab-demo.webm" type="video/webm" />
       <source src="http://rcg.group.shef.ac.uk/tutorial_videos/mobaxterm-login-matlab-demo.mp4" type="video/mp4" />
   </video>


What Next?
``````````

Now you have connected to a cluster, you can look at how to submit jobs with :ref:`sge-queue` or look at :ref:`sharc-software` and :ref:`iceberg-software`
