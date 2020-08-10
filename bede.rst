.. _bede:

Bede (Tier 2 GPU cluster)
=========================

Bede is a new `EPSRC-funded <https://gow.epsrc.ukri.org/NGBOViewGrant.aspx?GrantRef=EP/T022167/1>`__ 'Tier 2' (regional) HPC cluster.  
It is currently being configured and tested
and **should be available for use by researchers in the summer of 2020**.
The system will be available for use by researchers from `N8 Research Partnership`_ institutions
(Durham, Lancaster, Leeds, Liverpool, Manchester, Newcastle, Sheffield and York).

This system will be particularly well suited to supporting:
 
- **Jobs that require multiple GPUs and possibly multiple nodes**
- **Jobs that require much movement of data between CPU and GPU memory**

NB the system was previously known as NICE-19.

Hardware, OS and scheduler
--------------------------

* Main GPU nodes (32x) - each (`IBM AC922`_) node has

  * 2x IBM POWER9_ CPUs (and two NUMA nodes), with
  * 2x `NVIDIA V100`_ GPUs per CPU
  * Each CPU is connected to its two GPUs via high-bandwidth, low-latency NVLink interconnects
    (helps if you need to move lots of data to/from GPU memory)

* Inference GPU nodes (6x `IBM IC922`_)
   
  * 4x nodes have `NVIDIA T4`_ inference GPUs 
  * 2x nodes have FPGAs instead (`Bitware 250-SoC`_, which is a Xilinx Zinq Ultrascale+ FPGA that has ARM on the package too)

* Networking

  * 100 Gb EDR Infiniband
    (high bandwith and low latency to support multi-node jobs)
  * 10 Gb Ethernet as a backup

* Storage: Lustre parallel file system (available over Infiniband and Ethernet network interfaces)
* Scheduler: Slurm

Research software
-----------------

The set of research software available on the cluster is yet to be finalised by is likely to include the following (subject to change):

* IBM Watson Machine Learning Community Edition

  * IBM Distributed Deep Learning (DDL)

    * Efficiently scale popular machine learning frameworks over multiple CPUs/GPUs/nodes
    * Works with TensorFlow, IBMCaffe and Pytorch

  * IBM Large Model Support

    * Work with models too large to fit into the memory of a single GPU by transparently moving data from CPU to GPU memory as required.
    * Works with BVLC Caffe, IBMCaffe, TensorFlow, TensorFlow-Keras, PyTorch

  * IBM Snap ML

    * A library for training generalized linear models
    * Supports GPU acceleration, distributed training and sparse data structures
    * Can integrate with Scikit-Learn and Apache Spark

  * IBM provide this software inc their custom versions of PyTorch, TensorFlow etc via conda channels

* Profiling and debugging

  * Standard GNU toolkit via the IBM Advanced Toolchain for Linux

    * Provides IBM-optimised GNU compilers, BLAS/LAPACK, glibc, gdb, valgrind, itrace, Boost, Python, Go and more

  * NVIDIA tools

    * ``nvprof``
    * ``nsight-systems`` and ``nsight-compute``
    * ``cuda-gdb``

Acknowledgement for presentations and papers
--------------------------------------------

If you are hoping to use Bede in your research in the coming months then you will need to acknowledge the facility in presentations and papers using the following text:

   "We acknowledge the N8 Centre for Computationally Intensive Research (N8 CIR) for computational resources established through EPSRC (EP/T022167/1)"

A high-resolution version of the N8 CIR logo for use in talks or on posters can be found :ref:`here <N8 CIR logo>`.


.. _Bitware 250-SoC: https://www.bittware.com/fpga/250-soc/
.. _IBM AC922: https://www.ibm.com/uk-en/marketplace/power-systems-ac922
.. _IBM IC922: https://www.ibm.com/uk-en/marketplace/power-system-ic922
.. _N8 CIR logo: https://n8cir.org.uk/about/n8-cir-logo/
.. _N8 Research Partnership: https://www.n8research.org.uk/
.. _NVIDIA T4: https://www.nvidia.com/en-gb/data-center/tesla-t4/
.. _NVIDIA V100: https://www.nvidia.com/en-us/data-center/v100/
.. _POWER9: https://www.ibm.com/uk-en/it-infrastructure/power/power9
