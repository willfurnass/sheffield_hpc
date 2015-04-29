.. _cluster-specs:

Cluster Specifications
======================


Summary of iceberg hardware specs.
``````````````````````````````````

+ Total CPUs: 3440 cores
+ Total GPUs: 16 units
+ Total Memory: 31.8 TBytes
+ Permanent Filestore: 45 TBytes
+ Temporary Filestore: 260 TBytes
+ Physical size: 8 Racks
+ Maximum Power Consumption: 83.7 KWatts
+ All nodes are connected via fast infiniband.

For reliability, there are two iceberg head-nodes
'for loging in' configured to take over from each other
in case of hardware failure.


Worker Nodes CPU Specifications
```````````````````````````````

Intel Ivybridge based nodes


+ 92 nodes, each with 16 cores and 64 GB of total memory (i.e. 4 GB
  per core).
+ 4 nodes, each with 16 cores and 256 GB of total memory (i.e. 16GB
  per core).
+ Each node uses2 of Intel E5 2650V2 8-core processors (hence 2*8=16
  cores).
+ Scratch space on local disk of on each node is 400 GB


Intel Westmere based nodes


+ 103 nodes, each with 12 cores and 24 GB of total memory ( i.e. 2 GB
  per core )
+ 4 nodes with 12 cores and 48 GB of total memory ( i.e. 4GB per core
  )
+ Each node uses 2 of Intel X5650 6-core processors ( hence 2*6=12
  cores )


GPU Units Specifications
````````````````````````

8 Nvidia Tesla Kepler K40Ms GPU units


+ Each GPU unit contains 2880 thread processor cores
+ Each GPU unit has 12GB of GDR memory. Hence total GPU memory is
  8*12=96 GB
+ Each GPU unit is capable of about 4.3TFlop of single precision
  floating point performance, or 1.4TFlops at double precision.
+ Hencemaximum GPU processing power of is 11.2 TFlops in total.


8 Nvidia Tesla Fermi M2070s GPU units


+ Each GPU unit contains 448 thread processor cores
+ Each GPU unit contains 6GB of GDR memory. Hence total GPU memory is
  8*6=48 GB
+ Each GPU unit is capable of about 1TFlop of single precision
  floating point performance, or 0.5TFlops at double precision.
+ Hence maximum GPU processing power is 8 TFlops in total.


Filestore
`````````

+ 45 TBytes NFS mounted filestoreproviding users with storage on /home
  and /data areas
+ 260 TBytes Infiniband connected parallel filestore providing storage
  on /fastdata area

Filestore Allocations
`````````````````````
By default users get;


+ 5 GBytes of storage on their /home area
+ 50 GBytes of storage on their /data area.
+ Currently we set no limits on the /fastdata area but the files that
  have not been modified for 3 months will get deleted.


From time to time we shall review our data storage and allocation
policies depending on usage and inform our users.

It is strongly recommended that anyone wishing to use the /fastdata
area creates a subdirectory with the same name as their username for
storing their data.

The /fastdata area is faster to access from the new (INTEL-based)
nodes where the new infiniband connections are in use.


Filestore recovery and backup policies
``````````````````````````````````````

If you do loose files by accidental deletion, over-writing etc. and
you wish us to recover them, do let us know as soon as possible to
increase the chances of successful recovery.


+ Users' /home areas are fully backed up to allow recovery of lost
  data.
+ /data and /fastdata areas are not backed up, however ...
+ Due to mirroring it will usually be possible to recover lost or
  deleted files from the /data areas, providing we are informed quickly
  after such an incident.
+ It is not possible to recover lost and/or deleted files from the
  /fastdata areas.

Software and Operating System
`````````````````````````````

Users normally log into a head node and then use one (or more) of the
worker nodes to run their jobs on. Scheduling of users' jobs on the
worker nodes are managed by the 'Sun Grid Engine' software. Jobs can
be run interactively ( qsh ) or submitted as batch jobs ( qsub ).


+ The operating system is 64-bit Scientific Linux "which is Redhat
  based" on all nodes
+ The Sun Grid Engine for batch and interactive job scheduling
+ Many Applications, Compilers, Libraries and Parallel Processing
  Tools. See Section on Software.

