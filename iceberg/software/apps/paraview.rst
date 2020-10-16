.. Warning:: 
    Iceberg reaches end-of-life on **30th November 2020.**
    If you are running jobs on Iceberg then you need to take urgent action to ensure that your jobs/scripts will run on ShARC or Bessemer. 
 
    If you have never used ShARC or Bessemer then now is the time to test your scripts.
    Not all software on Iceberg is available on ShARC/Bessemer. 

Paraview
========

.. sidebar:: Paraview
   
   :Version: 4.3
   :Support Level: extras
   :Dependancies: openmpi (1.8)
   :URL: http://paraview.org/ 
   :Documentation: http://www.paraview.org/documentation/ 

Paraview is a parallel visualisation tool.

Using Paraview on iceberg
-------------------------

This guide describes how to use `Paraview <http://www.paraview.org/>`_ from iceberg.
Paraview is a parallel visualisation tool designed for use on clusters like iceberg.
Because iceberg is designed primarily as a headless computational cluster Paraview
has not been installed on iceberg in such a way that you load the GUI remotely on iceberg [1]_.
Paraview therefore runs on iceberg in client + server mode, the server running 
on iceberg, and the client on your local machine.  You need to: 

#. **Configure and start the client on your machine** first then
#. **then start the server on Iceberg**, which will then connect back to the client.

Configure the Client
^^^^^^^^^^^^^^^^^^^^

#. Download and install `Paraview 4.3 <http://www.paraview.org/download/>`_ [2]_ on your own machine.
#. Start the ``paraview`` program on your own machine.
#. Click **File** -> **Connect** -> **Add Server**
#. **Name**: ``my_iceberg_paraview_config``
#. **Server Type**: **Client / Server (reverse connection)** (we want Iceberg to connect to our local machine)
#. **Port**: ``11111`` (the default)
#. Click **Configure** to go to the next screen
#. **Startup Type**: **Manual** (the default)
#. Click **Save**
#. Back at the original **Choose Server Configuration** menu, click **Connect** to start listening for a connection from Iceberg

.. note:: 
    This method requires that Iceberg can connect to port 11111 on your machine.  
    You **may need to modify your machine's firewall** to permit such connections.  
    On Linux machines using the `UFW <https://wiki.archlinux.org/index.php/Uncomplicated_Firewall>`_ firewall you can allow connections
    to port 11111 on your machine from other machines on the University network (including Iceberg) and from no other machines using: ::

            sudo ufw allow from 143.167.0.0/16 to any port 11111

Start the Server
^^^^^^^^^^^^^^^^

After configuring the client:

#. Log in to iceberg **from the client machine** via ssh [3]_ 
#.  Run ``qsub-paraview`` from this login node (**not** from a worker node)

This will submit a job to the scheduler queue for 16 processes with 4GB of RAM each.
This is designed to be used for large visualisation tasks, smaller jobs can be 
requested by specifying standard ``qsub`` commands to ``qsub-paraview``.  For example, 
to request just one process: ::

        qsub-paraview -pe openmpi-ib 1

Assuming you still have the client listening for connections, once the paraview
job starts in the queue it should connect to your client and you should be able 
to start accessing data stored on iceberg and rendering images.

.. [1] It is not possible to install the latest version of the paraview GUI on  
   iceberg due to the Qt version shipped with Scientific Linux 5.
.. [2] The client and server versions have to match.
.. [3] Connecting to Paraview via the automatic method descibed here is not 
   supported on the MyApps portal.

A Note on Performance
---------------------

When you run Paraview locally it will use the graphics hardware of your local 
machine for rendering. This is using hardware in your computer to create and 
display these images. On iceberg there is no such hardware, it is simulated via
software. Therefore for small datasets you will probably find you are paying a 
performance penalty for using Paraview on iceberg.

The advantage however, is that the renderer is on the same cluster as your data,
so no data transfer is needed, and you have access to very large amounts of 
memory for visualising very large datasets.


Manually Starting the Server
----------------------------
The ``qsub-paraview`` command is a wrapper that automatically detects the client
IP address from the SSH connection and submits the job.
It is possible to customise this behavior by copying and modifying this script.
This for instance would allow you to start the paraview server via MyApps or 
from a different computer to the one with the client installed.
The script used by ``qsub-paraview`` also serves as a good example script and 
can be copied into your home directory by running ``cp /usr/local/bin/pvserver_submit.sh ~/``.
This script can then be qsubmitted as normal by ``qsub``.
The client IP address can be added manually by replacing ``echo $SSH_CLIENT | awk '{ print $1}'``
with the IP address.
More information on Paraview client/server can be found 
`here <http://www.paraview.org/Wiki/Setting_up_a_ParaView_Server#Running_the_Server>`_.


Installation
------------

Custom build scripts are availible in ``/usr/local/extras/paraview/build_scripts``
which can be used to recompile.
