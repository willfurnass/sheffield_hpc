.. _GPUIntro:

GPU Computing
=============

Graphics Processing Units (GPUs) were, as the name suggests, originally designed for the efficient processing of graphics.
Over time, they were developed into systems that were capable of performing general purpose computing
which is why some people refer to modern GPUs as GP-GPUs (General Purpose Graphical Processing Units).

Graphics processing typically involves relatively simple computations
that need to applied to millions of on-screen pixels in parallel.
As such, GPUs tend to be very quick and efficient at computing certain types of parallel workloads.

The `GPUComputing@Sheffield website
<http://gpucomputing.shef.ac.uk/>`_ aims to facilitate the use of GPU computing within University of Sheffield research by
providing resources for training and purchasing of equipment as well as providing a network of GPU users and research projects within the University.

.. _GPUCommunity:

GPU Community and NVIDIA Research Centre
----------------------------------------
The University of Sheffield has been officially affiliated with `NVIDIA <https://research.nvidia.com/>`_ since 2011
as an `NVIDIA CUDA Research Centre <https://developer.nvidia.com/academia/centers/university-sheffield-cuda-research-center>`_.
As such, NVIDIA offer us some benefits as a research institution including
discounts on hardware,
technical liaisons,
online training and
priority seed hardware for new GPU architectures.
For first access to hardware, training and details of upcoming events, discussions and help
please join the `GPUComputing google group <https://groups.google.com/a/sheffield.ac.uk/forum/#!forum/gpucomputing>`_.

.. _GPUAccess:

Requesting access to GPU facilities
-----------------------------------

**Public GPU nodes are available to in Bessemer, ShARC and Iceberg. These can be be used without acquiring extra permission.**

Research groups also have an option to purchase and add nodes to the cluster to be managed by IT Services.
For these nodes (e.g. :ref:`dcs_gpu_nodes_bessemer`),
permission from the owner(s) of those nodes is required for access.

The node owner(s) always have highest priority on the job queue but
as a condition for the management of additional nodes on the cluster,
the nodes are expected to be used as a resource for running short jobs during idle time.
If you would like more information about having IT Services add and manage custom nodes,
please contact ``research-it@sheffield.ac.uk``.

Using GPUs on the different clusters
------------------------------------

Hardware, software and ways to access GPUs differ between the clusters.
Instructions specific to each cluster can be found below:

* :ref:`GPUComputing_bessemer`
* :ref:`GPUComputing_sharc`
* :ref:`GPUComputing_iceberg`
