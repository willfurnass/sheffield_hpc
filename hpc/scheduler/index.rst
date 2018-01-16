.. _scheduler:

Using interactive sessions / submitting batch jobs
==================================================

.. toctree::
    :hidden:
    :maxdepth: 1
    :glob:

    ./*

Getting started 
---------------

See our guide to :ref:`sge-queue`.

Submitting and controlling jobs programmatically
------------------------------------------------

See our guide to the :ref:`drmaa` API.

Reference information
---------------------

Commands that allow you to interact with the scheduler:

* :ref:`qhost` - Show's the status of Sun Grid Engine hosts.
* :ref:`qrsh` - Requests an interactive session on a worker node. No support for graphical applications.
* :ref:`qrshx` - Requests an interactive session on a worker node. Supports graphical applications. Superior to :ref:`qsh` in most cases.
* :ref:`qsh` - Requests an interactive session on a worker node. Supports graphical applications.
* :ref:`qstat` - Displays the status of jobs and queues.
* :ref:`qsub` - Submits a batch job to the system.
* :ref:`qtop` - Provides a summary of all processes running on the cluster for a given user
