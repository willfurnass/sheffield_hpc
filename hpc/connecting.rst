.. _connecting:

Connecting to a cluster using SSH
=================================

The most versatile way to **run commands and submit jobs** on one of the clusters is to 
use a mechanism called `SSH <https://en.wikipedia.org/wiki/Secure_Shell>`__, 
which is a common way of remotely logging in to computers 
running the Linux operating system.  

To connect to another machine using SSH you need to 
have a SSH *client* program installed on your machine.  
macOS and Linux come with a command-line (text-only) SSH client pre-installed.  
On Windows there are various graphical SSH clients you can use, 
including *MobaXTerm*.

**Whether/how you can connect** to a University cluster using SSH (or the related protocols SCP and SFTP) **depends on**:

* Where you are connecting from:

  * on campus;
  * off campus *after* `establishing a VPN connection <https://www.sheffield.ac.uk/it-services/vpn>`_;
  * off campus without a VPN connection.

* Whether `Multifactor Authentication (MFA) <https://sites.google.com/sheffield.ac.uk/mfa/home>`__  has been enabled on your University account.
  If MFA has been enabled then when logging in with SSH you will be prompted to enter a one-time code or send a push notication to your MFA device 
  after entering your username and password.  If you do not have MFA enabled on your account then you will not be able to login from off campus without using VPN.

* Whether you want to use password-based authentication or 'public-key'-based authentication.

**Authentication requirements per cluster**:

+----------+------------------------+---------------------------------------------------------------------------------------------------+
| Cluster  | From campus or via VPN | From off campus and without a VPN connection                                                      |
+==========+========================+===================================================================================================+
| Bessemer | Password or public key | Not permitted                                                                                     |
+----------+------------------------+---------------------------------------------------------------------------------------------------+
| ShARC    | Password or public key | Not permitted                                                                                     |
+----------+------------------------+---------------------------------------------------------------------------------------------------+

SSH client software on Windows
------------------------------

Download and install the *Installer edition* of `MobaXterm <https://mobaxterm.mobatek.net/download-home-edition.html>`_.

After starting MobaXterm you should see something like this:

.. image:: /images/mobaxterm-welcome.png
   :width: 50%
   :align: center

Click *Start local terminal* and if you see something like the following then please continue to :ref:`ssh`.

.. image:: /images/mobaxterm-terminal.png
   :width: 50%
   :align: center

Running commands from a terminal (from the command-line) may initially be
unfamiliar to Windows users but this is the recommended approach for
running commands on Bessemer or ShARC as 
it is the idiomatic way of interfacing with the Linux clusters.

SSH client software on Mac OS/X and Linux
-----------------------------------------

Linux and macOS (OS X) both typically come with a command-line SSH client pre-installed.

If you are using macOS and want to be able to run graphical applications on the clusters then 
you need to install the latest version of the `XQuartz <https://www.xquartz.org/>`_ *X Windows server*.

Open a terminal (e.g. *Gnome Terminal* on Linux or *Terminal* on macOS) and then go to :ref:`ssh`.

.. _ssh:

Establishing a SSH connection
-----------------------------

Once you have a terminal open run the following command to 
log in to a cluster: ::

    ssh -X $USER@$CLUSTER_NAME.shef.ac.uk

Here you need to:

* replace ``$USER`` with your IT Services username (e.g. ``te1st``)
* replace ``$CLUSTER_NAME`` with ``bessemer`` or ``sharc``.

.. note::

    **macOS users**: if this fails then:

    * Check that your `XQuartz <https://www.xquartz.org/>`_ is up to date then try again *or*
    * Try again with ``-Y`` instead of ``-X``

This should give you a prompt resembling the one below: ::

    [te1st@sharc-login1 ~]$

At this prompt type: ::

    qsh

Like this: ::

    [te1st@sharc-login2 ~]$ qsh
    Your job 135355 ("INTERACTIVE") has been submitted
    waiting for interactive job to be scheduled ....
    Your interactive job 135355 has been successfully scheduled.

Which will pop up another terminal window, which supports graphical applications.

.. note::

    When you login to a cluster you reach one of two login nodes. 
    You **should not** run applications on the login nodes.
    Running ``qsh`` gives you an interactive terminal 
    on one of the many worker nodes in the cluster.

    If you only need terminal-based (command-line only) applications 
    you can run the ``qrsh`` command,
    which will give you a shell on a worker node, 
    but without graphical application (X server) support.


.. raw:: html

   <p>
    This video shows the connection process using MobaXterm, and then connection
    and running MATLAB from a <code>qsh</code> terminal.
   </p>

   <video style="margin-left: auto; margin-right:auto; display: block;" width=70% controls>
       <source src="https://rcg.group.shef.ac.uk/tutorial_videos/mobaxterm-login-matlab-demo.webm" type="video/webm" />
       <source src="https://rcg.group.shef.ac.uk/tutorial_videos/mobaxterm-login-matlab-demo.mp4" type="video/mp4" />
   </video>

HPC SSH Gateway
---------------

Direct SSH access to the HPC clusters from off campus is not possible without the use of VPN. However 
if you are unable to use VPN we also provide an SSH gateway service to allow off-site SSH access to our HPC clusters.

.. note::

  Use of a `VPN connection <https://www.sheffield.ac.uk/it-services/vpn>`_ is the recommended method to use for off-site SSH access to the HPC clusters.


.. note::
  * Access to the HPC SSH gateway service requires that you have an existing :ref:`HPC account <accounts>`. 
  * You must additionally request access to the HPC SSH gateway by emailing `it-servicedesk@sheffield.ac.uk <it-servicedesk@sheffield.ac.uk>`_ including a justification for your request.

The SSH gateway servers are configured to be SSH jump hosts only.  They do not have direct access to HPC filestore or Research Shared areas, and you cannot 
run an interactive SSH terminal session directly on the gateway servers.  Additionally the HPC gateway servers only allow access to the HPC clusters - you cannot access any 
other IT Services or departmental servers using this gateway.


* Access a HPC cluster via SSH: ::

    ssh -J [username]@hpcgw.shef.ac.uk [username]@sharc.shef.ac.uk

* Transfer a file using SCP: ::

    scp -J [username]@hpcgw.shef.ac.uk [source path] [destination path]

* Transfer files using Rsync: ::

    rsync -av -e 'ssh -J [username]@hpcgw.shef.ac.uk' [source path] [destination path]


* Using WinSCP: ::

    New Session -> Advanced -> Connection -> Tunnel
    Select 'Connect through SSH tunnel'
    Hostname: 'hpcgw.shef.ac.uk'
    Port number: '22'

.. image:: /images/SSHgatewayWinSCP.png
   :width: 75%
   :align: center

* Configure MobaXterm: :: 

    Edit 'Session Settings':
    Set 'SSH Use 2-factor authentication for SSH gateways'

.. image:: /images/SSHgatewayMobaXtermSettings.png
   :width: 75%
   :align: center

* Create a new session using MobaXterm: ::

    Select 'Network settings' tab within SSH Session settings
    Select 'Connect through SSH gateway (jump host)
    Gateway SSH server: 'hpcgw.shef.ac.uk'
    Port: '22'

.. image:: /images/SSHgatewayMobaXtermSession.png
   :width: 75%
   :align: center

* When prompted to enter your Duo two-factor code either input a 6 digit code from your Duo device or enter '1' for a push notification to be sent to your device.

What Next?
----------

Now you have connected to a cluster, 
you can look at how to submit jobs with :ref:`submit-queue` or 
look at the software installed on 
:ref:`Bessemer <bessemer-software>`,
:ref:`ShARC <sharc-software>` and 
