.. _mfa:

Tips for working with MFA
=========================

Here are some techniques you can use to reduce the number of times that you need to reauthenticate to our HPC systems.

.. include:: /referenceinfo/imports/staying_connected.rst


TMUX/screen
-----------

`TMUX <https://github.com/tmux/tmux/wiki>`_ and `screen <https://www.gnu.org/software/screen/manual/screen.html>`_
are available on the HPC login nodes and can be used to run multiple shell sessions within a single SSH session.

.. note::

    A *shell* is the command-line interface where you type commands.
    Normally, when you SSH into the cluster you only get one shell session.

    Tools like TMUX and screen let you:

    - Detach from and reattach to sessions so your processes keep running even if you disconnect.
    - Run several virtual terminals concurrently.

For examples of using TMUX to manage multiple sessions see the following RSE blog post: 
`tmux: remote terminal management and multiplexing <https://rse.shef.ac.uk/blog/tmux-intro/>`_

For a quick reference, see `TMUX Cheat Sheet <https://tmuxcheatsheet.com/>`_.

.. tip::

   .. code-block:: bash

        tmux new-session -A -s main

   This command attaches to a session named "main" if it exists on the current node, or creates it if not.

   It's recommended to run this command as your first command when logging in via SSH,
   so you immediately attach to your persistent session.


Transferring files
------------------

If you need to transfer data from your local PC to a research shared directory you can directly access the data from
your local PC without using MFA, instead of transferring the files via the HPC.  

For more info on how to do this see our
`Research Storage documentation <https://www.sheffield.ac.uk/it-services/research-storage/using-research-storage>`_ .
