.. _stanage-connecting-ssh:

.. toctree::
    :hidden:
    :maxdepth: 1
    :glob:

Connecting with SSH
===================

.. Hint::

    Usernames to connect with all HPC services will be the same as those you use to login to MUSE :underline-bold:`not` the prefix on your email address.


Once you have a terminal open run the following command to
log in to a cluster: ::

    ssh -X $USER@stanage-login1.shef.ac.uk

Launch an interactive session on Stanage using the following command: ::

    srun --mem=XXG --pty bash -i

where XX is the memory request for the session (default is 2Gb).



