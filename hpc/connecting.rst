.. _connecting:

Connecting to a cluster using myApps (web browser)
==================================================

In order to access ShARC and Bessemer you must set up a 'VPN and MFA <https://www.sheffield.ac.uk/it-services/vpn>`_. Also see section **Whether/how you can connect** below. 

The web browser method of access to ShARC and Bessemer is provided by the 'Secure Global Desktop Client <https://myapps.shef.ac.uk/installic_en.html>`_.

This method works well for most browsers on all the common computing platforms (Linux, Windows, Mac), however we recommend Internet Explorer or Firefox.

To login to ShARC or Bessemer click the following link:

'Connect via myAPPs Portal <https://myapps.shef.ac.uk/sgd/index.jsp?langSelected=en>`_

Logging in for the first time:

When you login for the first time, select Client Options on the myApps Portal page (bottom right).

Select Install the Client for the top option To use the installed Oracle Secure Global Desktop Client.

This will take you to the client download/install page.

Download and install the client for your system (Windows, Mac, or Linux).

Enter your username and password on the myApps Portal login page

Once you have managed to login, you will see a window with Applications on the Left Hand Panel.

there are icons for ShARC Applications & Bessemer Applications respectively.

For each of these you can select a HPC interactive job or a HPC terminal (where HPC is ShARC & Bessemer).
The interactive job is equivalent to a qsh/qrshx session on a worker node
The terminal is equivalent to a login node session from which you can use qsh/qrshx, qrsh, qsh-vis


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
  If MFA has been enabled then when logging in with SSH you will be prompted to enter a one-time code or send a push notification to your MFA device
  after entering your username and password.  If you do not have MFA enabled on your account then you will not be able to login from off campus without using VPN.

* Whether you want to use password-based authentication or 'public-key'-based authentication.

.. note::

    The `University Connect for China (UCC) <https://www.sheffield.ac.uk/it-services/university-connect-china>`_ is not the same service as the SSL VPN service and will not grant access to the HPC clusters.
    Users of the UCC must disconnect the UCC and connect to the SSL VPN in order to connect to the HPC clusters.

**Authentication requirements per cluster**:

+----------+------------------------+---------------------------------------------------------------------------------------------------+
| Cluster  | From campus or via VPN | From off campus and without a VPN connection                                                      |
+==========+========================+===================================================================================================+
| Bessemer | Password or public key | Not permitted (unless using the :ref:`HPC SSH gateway service <hpcgw_summary>`)                   |
+----------+------------------------+---------------------------------------------------------------------------------------------------+
| ShARC    | Password or public key | Not permitted (unless using the :ref:`HPC SSH gateway service <hpcgw_summary>`)                   |
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

---------

.. _hpcgw_summary:

What if I cannot use the VPN or I need a persistent long term connection
---------------------------------------------------------------------------

Direct SSH access to the HPC clusters from off campus is not possible without the use of VPN. However
if you are unable to use VPN we also provide an SSH gateway service to allow off-site SSH access to our HPC clusters.

.. note::
  * Access to the HPC SSH gateway service requires that you have an existing :ref:`HPC account <accounts>`.
  * You must additionally request access to the HPC SSH gateway by emailing `it-servicedesk@sheffield.ac.uk <it-servicedesk@sheffield.ac.uk>`_ including a justification for your request.
  * If the cluster access can be handled via the usage of the SSL VPN without undue effort, your request will not be granted.

For more information see :ref:`HPC Gateway Service Details <hpcgateway>`.



What Next?
----------

Now you have connected to a cluster,
you can look at how to submit jobs with :ref:`submit-queue` or
look at the software installed on
:ref:`Bessemer <bessemer-software>` and
:ref:`ShARC <sharc-software>`
