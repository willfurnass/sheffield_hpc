.. _stanage-file-transfers:

.. toctree::
    :hidden:
    :maxdepth: 1
    :glob:

Transfering Files
=================

.. warning::

  As with all connections to the clusters, if you are not using a wired ethernet connection in a 
  University campus building then you will need to `turn on the VPN <https://www.sheffield.ac.uk/it-services/vpn>`_.

To transfer files to/from the clusters you can:

* Use a program that supports one or both of the **SCP** and **SFTP** protocols to copy/move files to/from your own machine 
  or from a remote machine to the cluster.
* Use a program like ``curl`` or ``wget`` to download files directly to the clusters.

---------


Transfers with SCP/SFTP
-----------------------

Secure copy protocol (SCP) is a protocol for securely transferring computer files between a local host and a 
remote host or between two remote hosts. It is based on the Secure Shell (SSH) protocol and the acronym typically 
refers to both the protocol and the command itself.

Secure File Transfer Protocol (SFTP) is also a file transfer protocol. It is based on the 
FTP protocol with included SSH security components.

.. hint::

  If you need to move large files (e.g. larger than a gigabyte) from one remote machine to the cluster you 
  should SSH in to the computer hosting the files and use scp or rsync to transfer over to the other directly as this will 
  usually be quicker and more reliable.

.. raw:: html

    <hr class="hr-mid-section-separator">

Using SCP in the terminal
^^^^^^^^^^^^^^^^^^^^^^^^^

If your local machine has a terminal and the ``scp``  (“secure copy”) command is available 
you can use it to make transfers of files or folders.

Where below substitute **$USER** with your cluster username. 
You should be prompted for your Duo MFA credentials after entering your password. Request a push notification or enter your passcode.

To upload, you transfer from your local machine to the remote cluster:

.. code-block:: shell

  scp /home/user/file.txt $USER@stanage-login1.shef.ac.uk:/users/$USER/

To download, you transfer from the remote cluster to your local machine:

.. code-block:: shell

  scp $USER@stanage-login1.shef.ac.uk:/home/$USER/file.txt /home/user/

To copy a whole directory, we add the ``-r`` flag, for “recursive”

.. code-block:: shell

  scp -r $USER@stanage-login1.shef.ac.uk:/home/$USER/my_results /home/user/


.. raw:: html

    <hr class="hr-mid-section-separator">

Using Filezilla
^^^^^^^^^^^^^^^^^^^^

FileZilla is a cross-platform client available for Windows, MacOS and Linux for downloading 
and uploading files to and from a remote computer.

Download and install the FileZilla **client** from https://filezilla-project.org. After installing and opening the program, 
there is a window with a file browser of your local system on the left hand side of the screen
and when you connected to a cluster, your cluster files will appear on the right hand side.

To connect to the cluster, we’ll just need make a **new site** and enter our credentials in the **General** tab:

.. caution::

  By default Filezilla will save profiles in plaintext on your machine. You must ensure you use a master password to 
  encrypt these credentials by changing the settings 
  `as shown in these instructions <https://filezillapro.com/docs/v3/advanced/master-password/>`_.

* **Host**: sftp://stanage-login1.shef.ac.uk
* **User**: Your cluster username
* **Password**: Your cluster password (leave blank and fill this interactively if on a shared machine.)
* **Port**: (leave blank to use the default port)
* **Protocol**: sftp
* **Logon Type**: Interactive

In the **transfer settings** tab limit the number of simultaneous connections to 1.

Save these details as a profile and then connect. You should be prompted for your Duo MFA credentials. 
Request a push notification or enter your passcode.  You will now see your remote files appear on the 
right hand side of the screen. This process can be repeated to save a profile for each cluster.

You can drag-and-drop files between the left (local) and right (remote) sides of the screen to transfer files.

.. raw:: html

    <hr class="hr-mid-section-separator">

Using rsync
^^^^^^^^^^^^^^^^^^^^

As you become more familiar with transferring files, you may find that the ``scp`` is limited. The ``rsync`` utility provides 
advanced features for file transfer and is typically faster compared to both ``scp`` and ``sftp``. It is a utility for 
efficiently transferring and synchronizing files between storage locations including networked computers by comparing the 
modification times and sizes of files. The utility is particularly useful as it can also resume failed or partial file 
transfers by using the ``--append-verify`` flag.

Many users find ``rsync`` is especially useful for transferring large and/or many files as well as creating synced 
backup folders.

.. caution::

  It is easy to make mistakes with ``rsync`` and accidentally transfer files to the wrong location, sync in the wrong 
  direction or otherwise accidentally overwrite files. To help you avoid this, you can first use the ``--dry-run`` flag for 
  ``rsync`` to show you the changes it will make for a given command.

The ``rsync`` syntax is very similar to ``scp``. To transfer to another computer with commonly used options, 
where below substitute **$CLUSTER_NAME** with bessemer or sharc and **$USER** with your cluster username.
You should be prompted for your Duo MFA credentials after entering your password. Request a push notification or 
enter your passcode:

.. code-block:: shell

  rsync -avzP /home/user/file.iso $USER@stanage-login1.shef.ac.uk:/home/$USER/

The ``a`` (archive) option preserves file timestamps and permissions among other things; 
the ``v`` (verbose) option gives verbose output to help monitor the transfer; 
the ``z`` (compression) option compresses the file during transit to reduce size and transfer time; 
and the ``P`` (partial/progress) option preserves partially transferred files in case of an interruption 
and also displays the progress of the transfer.

To recursively copy a directory, we can use the same options:

.. code-block:: shell

  rsync -avzP /home/user/isos/ $USER@stanage-login1.shef.ac.uk:/home/$USER/

This will copy the local directory and its contents under the specified directory on the remote system. 
If the trailing slash is omitted on the destination path, a new directory corresponding to the transferred 
directory (isos in the example) will not be created, and the contents of the source directory will be copied 
directly into the destination directory.

As before with ``scp``, to download from the cluster rather than upload simply reverse the source and destination:

.. code-block:: shell

  rsync -avzP $USER@stanage-login1.shef.ac.uk:/home/$USER/isos /home/user/ 

---------
