.. _IcebergCosts:

.. toctree::
   :maxdepth: 1
   :hidden:


Iceberg Resource Costs
======================
Normal usage of centrally provided HPC facilities 'i.e. iceberg' at the University of Sheffield is free. For research grants, as with other CICS services, it is included in the central service overhead. This means that users receive a fair share of HPC services alongside other users, who number about 900 in any one academic year.

Researchers requiring compute resources beyond the capabilities of  "free at the point of use" service have the following options:

* :ref:`Use the N8 Regional Computing Facilities <UsetheN8Facilities>`
* :ref:`Purchase HPC hardware <PurchaseHPCHardware>`
* :ref:`Reserve a portion of the local HPC facilities for a length of time <ReserveHPC>`

Use of the N8 Facilities
------------------------

.. _UsetheN8Facilities:

The University of Sheffield has a 1/8 share of the N8 HPC facility and manages access to it on behalf of researchers at The University of Sheffield. Details on how to access the facility are at http://n8hpc.org.uk/getting-started

We also provide some local N8 related information on our web pages. FIXME

Purchase of HPC Hardware
------------------------

.. _PurchaseHPCHardware:

Research groups can work with CICS to purchase hardware that will sit within the framework of iceberg while delivering dedicated HPC services for their research. This will give research teams the advantage of having access to dedicated resources while continuing to take advantage of the 'free' facilities as well. The University of Sheffield has a framework agreement for procuring such extra HPC hardware at favourable prices. This is prepared by CICS working with the Research Computing Advisory Group. Research groups who are interested should initially contact hpchub@sheffield.ac.uk to start a dialog.  
It is important to note that the 'framework agreement' ensures and in a way limits the choice of new hardware so as to be able to integrate it with the iceberg cluster.

Reservation of HPC Compute Resources
------------------------------------

.. _ReserveHPC:

Researchers can reserve a fraction of the central HPC facility for a length of time which they specify. This is paid for using a project research grant or other funds. See current HPC purchasing costs and procedures for more details.

Extra File Storage Costs
------------------------

CiCS also provides general purpose networked file storage at prices ranging from £260 per TB per annum to £900 per TB per annum depending on backup requirements (e.g if there is a requirement to backup to tape, or mirror data to a separate location) - for more information on this please see the `CICS information on filestore <http://www.shef.ac.uk/cics/filestore>`_.

Current HPC Purchasing and costs
--------------------------------
The current costs for dedicated HPC resources are shown in the tables below.  

To apply for dedicated resource, researchers should e-mail hpchub@sheffield.ac.uk and provide the following information:

* Project name
* Contact details for the project which will be paying for the resource
* The resource requirements of the project
* The start date of the reservation
* How long the reservation is for

Researchers should also let us know if they require reservations for licensed commercial packages and if dedicated support from CICS research computing staff is required. It will be helpful to specify finer details of the compute requirements, such as-  if you are running a parallel application, how many cores and how much memory is required. Is the application shared memory or distributed memory.

We will reply and give you details of the cost, how to pay and how to access your reserved allocation. As it is necessary to reserve resources the process can take a week (the HPC facility is well utilised) it can take longer if more specialised and/or a greater number of resources is requested. We will inform you as soon as the reservation is active and you will receive a warning of when the reservation is due to expire.The costs quoted below are based on the full economic cost and include staffing costs for maintaining the facility. Income generated using these resource allocation requests is reinvested in facility upgrades.

Costs per CORE per hour:

========== ====================== ==========================================================
Resource   Cost per CORE per hour Notes
---------- ---------------------- ----------------------------------------------------------
standard   0.98 pence             Single Intel E5-2650 v2 core with 4 GB per core of memory
---------- ---------------------- ----------------------------------------------------------
Big Memory 1.13 pence             Single Intel E5-2650 v2 core with 16 GB per core of memory
========== ====================== ==========================================================

Costs per NODE per hour:

========== ====================== ==========================================================
Resource   Cost per NODE per hour Notes
---------- ---------------------- ----------------------------------------------------------
Standard   15.8 pence             Dual Intel E5-2650 v2 8-core CPU with 4GB/core
---------- ---------------------- ----------------------------------------------------------
Big Memory 18.2 pence             Dual Intel E5-2650 v2 8-core CPU with 16GB/core
---------- ---------------------- ----------------------------------------------------------
GPU        38.4 pence             2 x NVIDIA Kepler K40M GPUs with 2xE5-2650-v2 CPUs
========== ====================== ==========================================================

