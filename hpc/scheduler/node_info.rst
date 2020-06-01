Cluster information
===================

You may occasionally want to see information on what nodes are available and on the used and free resources associated with those nodes.  
The commands here may help with that.

SGE
---

.. _qhost:

qhost
^^^^^

`qhost` is a SGE command that shows the state of worker nodes inc. available and used resources.

Run it using just the following from any login or worker node: ::

    qhost

By default this shows every node in the cluster: note that some of these nodes may be reserved/purchased for specific research groups and hence may not be available for general use.

To filter the ``qhost`` output or get additional information from it read up on its usage by running: ::

    man qhost

Slurm
-----

sinfo
^^^^^

Show info on all Slurm *partitions* (worker node groups): ::

   sinfo --long

Show per-node info for all Slurm partitions: ::

   sinfo --long --Node
