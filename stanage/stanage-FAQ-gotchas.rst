.. _stanage-FAQ-gotchas:

Stanage FAQ and Gotchas
=======================

This section contains important information about the usage of the Stanage cluster and any significant changes
from expected behaviour on previous University of Sheffield HPC clusters.

-----

Graphical applications
----------------------

At present, we only provide very limited support for running graphical sessions on Stanage. Please see :ref:`Graphical sessions on Stanage<flight-desktop>`

-----

Poor performance with OpenMPI 4.1.1
-----------------------------------

Current installations of OpenMPI 4.1.1 have poor performance over Omnipath which is still under investigation. 
These versions of OpenMPI will be in use for the ``foss-2021`` Easybuild toolchain and should be avoided.

-----

Shared research area directories
--------------------------------

Shared research areas are available on the Stanage cluster :underline-bold:`upon request`. For details on this please see our page on :ref:`shared_dir`.

-----

My tmux session is missing from the login node?  
-----------------------------------------------

The Stanage cluster uses two login nodes in a active-backup configuration. This means when you connect via ``stanage.shef.ac.uk`` you will actually connect to
one of the two login nodes, ``stanage-login1.shef.ac.uk`` or ``stanage-login2.shef.ac.uk`` depending on which is currently the active login node.

If you disconnect from the cluster and return at a later time, you may not connect back to the same login node and won't be able to access the session you left running. 
If you are going to make use of ``tmux`` sessions on Stanage, you should choose one of the two login nodes and consistently directly connect via it.

It is also possible that after a maintenance period or a login node reboot, that all running ``tmux`` sessions will have been terminated as part of this 
and you will need to recreate your sessions. If you think this might have happened then you can check if the given login node's uptime with the ``uptime`` command to 
see if there has been a recent reboot.


