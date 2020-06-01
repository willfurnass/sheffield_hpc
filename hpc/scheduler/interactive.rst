.. _sched_interactive:

Interactive sessions
====================

*On SGE and SLURM*

If you wish to use a cluster for interactive work,
such as running applications like MATLAB or ANSYS,
or compiling software,
you will need to request an *interactive session* on a worker node from the scheduler.

This typically must be done from a cluster *login node*.
See :ref:`getting-started` for information on connecting to the clusters.

.. contents::
   :local:

On SGE
------
When using a :ref:`cluster that runs the SGE job scheduler <sched>`
you start an interactive shell session on a worker node using one of three commands,
``qrsh``, ``qsh`` or ``qrshx``,
plus, optionally, additional parameters to request resources
(e.g. run time, CPU cores, memory, GPUs, access to private worker nodes)
for your session.
For example, to start an interactive shell session on a worker node
with 2 CPU cores and 4 GB RAM per CPU core
SSH to a worker node then
run something similar to the following to start your interactive session: ::

    qrshx -pe smp 2 -l rmem=4G

The three possible commands differ as follows:

.. list-table::
   :widths: 25 75

   * - ``qrsh``
     - No support for graphical applications. Standard SGE command.
   * - ``qsh``
     - Supports graphical applications. Standard SGE command.
   * - ``qrshx``
     - Supports graphical applications. Superior to ``qsh``. Unique to Sheffield's clusters.

On SLURM
--------
When using a :ref:`cluster that runs the SLURM job scheduler <sched>`
an interactive shell session on a worker node is started by running the `srun <https://slurm.schedmd.com/srun.html>`__ command
on a login node
using something like: ::

   srun OPTIONS --pty /bin/bash -i

where ``OPTIONS`` is an optional set of parameters for requesting resources for the session
(run time, CPU cores, memory, GPUs, access to private worker nodes). 
For example, to start an interactive session with 4 CPU cores and 4GB of RAM per CPU core: ::

   srun OPTIONS -c4 --mem-per-cpu=4G --pty /bin/bash -i


Customising your job request using parameters
---------------------------------------------
Both SGE and SLURM support a **large number of parameters** for
**customising the resource requests associated with each interactive session**
(and :ref:`batch job <sched_batch>`).
As shown above, multiple parameters can be specified when starting an interactive session.
See :ref:`sched_options` for information on the most common/valuable parameters.

If you don't explicitly request resources such as run time, CPU cores or memory then
**default values** will be used.
For interactive sessions the defaults are typically:

* Run time: 8 hours
* CPU cores: 1
* RAM: 2GB

When to submit batch jobs instead
---------------------------------

.. note::

    Long running jobs *should* use the batch submission system
    rather than requesting an interactive session for a long time.
    Doing this will lead to better cluster resource utilisation for all users.
