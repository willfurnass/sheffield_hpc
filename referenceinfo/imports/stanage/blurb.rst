Stanage, will be the University's newest HPC system. It will be well suited to a variety of research workloads, 
including large jobs or jobs requiring GPUs. It is intended to replace the ShARC HPC system, but is much larger.

Stanage will have a similar number of worker nodes (192) to the existing ShARC cluster but these 
will each have more CPU cores (64 per node for a total ~12032 cores!), with a standard memory size of 256 GB.

The standard nodes have 256GB of memory per node (4GB/core) with an additional 12 large memory nodes with 
1TB  (16GB/core), and a further 12 very large nodes with 2TB (32GB/core).


TO BE CONFIRMED >>>

    The system as a whole has a larger performant shared storage area (4 PB traditional lustre filesystem) and 
    a novel 200 TB very fast storage area (Lustre filesystem on NVMe storage) optimised for random read/writes and 
    small file operations for temporary file storage.

<<< TO BE CONFIRMED

The system will also eventually feature 18 GPU nodes each with 4x NVIDIA A100 GPUs which will support machine 
learning research and other workloads that can benefit from GPU acceleration. 

.. note::

   At present these GPUs are currently installed in the Bessemer cluster.

Stanage uses the SLURM scheduler which is also in use on the Bessemer cluster.
