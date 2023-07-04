.. _sharc_decommissioning:

ShARC decommissioning
=====================

**ShARC will be decommissioned on 30th November 2023, after which time users will no longer be able to access that cluster and any jobs running or 
queueing at that time will be cancelled.**

Once ShARC has been decommissioned and all data securely deleted the system is to be donated to the `Taras Shevchenko National University of Kyiv <http://www.univ.kiev.ua/en/>`_.

We recommend users move to the :ref:`Stanage<stanage>` cluster which is a new, much larger system than previous TUOS clusters with drastically larger numbers of CPU cores, numbers of cores per node, memory per node (especially with Stanage's 1TB and 2TB nodes) and a large number of GPUs.

Please get in touch with us if you have any concerns regarding this announcement via HPC drop-in sessions (every Friday) or `research-it@sheffield.ac.uk. <mailto:research-it@sheffield.ac.uk?subject=ShARC%20HPC%20decommissioning>`_


Key things you need to know
---------------------------

.. :: 

----

Data must be migrated
"""""""""""""""""""""

Any user data stored in **/home**, **/data** or **/fastdata** directories on ShARC that you want to retain will :underline-bold:`need` to be migrated elsewhere e.g.

* to equivalent areas on Bessemer or Stanage or

* to **/shared/GROUPNAME** areas, which can be made accessible from Bessemer or from Stanage login nodes if they aren't already.

For guidance on moving your data please see our page on :ref:`transferring files<transferring_files>`.

----

ShARC specific shared filestorage areas will no longer be accessible
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

Shared application installations within **/usr/local** on ShARC (including **/usr/local/packages** and **/usr/local/community**) will not be accessible after the decommissioning date.

Equivalent shared application installation areas will be made available on the Stanage cluster prior to the ShARC decommissioning.

----

Stanage is our recommended target for migration 
"""""""""""""""""""""""""""""""""""""""""""""""

Stanage is our newest HPC cluster and is the logical successor to ShARC; it is well suited to running larger parallel jobs spanning multiple nodes and benefits from a high-performance network interconnect between nodes.  

As Bessemer only permits single node execution and has approximately 18 months of life left, we recommend that ShARC users preferentially migrate workloads to Stanage instead.

----

User installed software may require recompiling
"""""""""""""""""""""""""""""""""""""""""""""""

If you have compiled and installed your own applications on ShARC you should recompile your applications on the newer clusters when migrating your workloads in order to ensure compatibility and potentially significant performance increases. 

If you currently use IT Services-managed applications on ShARC which are not currently available on Stanage or Bessemer then you can submit a request to `research-it@sheffield.ac.uk <mailto:research-it@sheffield.ac.uk?subject=HPC%20Software%20installation%20request>`_.


