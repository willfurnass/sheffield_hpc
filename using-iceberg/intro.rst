.. _iceberg-tasks:

Types of Tasks on iceberg
=========================

The iceberg computing cluster is a very powerful and flexible service. 
It can be used for running many different types of computationally expensive 
tasks. The best way for you to use iceberg, and the tools you will need to learn
depend on what you want to do.


* Does your task take less than 8 hours to run?

If you have a 'short' running task, which takes less than 8 hours, either 
a graphical application such as MATLAB or Ansys or a script, you can 
use an interactive session. See :ref:`sge-interactive`.

* Does your task take longer than 8 hours, or do you have lots of short tasks?

If you need more than 8 hours of compute time, or you need to run lots of 
small tasks you will need to learn the Grid Engine submission system.
See :ref:`sge-intro`.

* Does your graphical or interactive taks require multiple threads or lots of memory?

You will need to specify Grid Engine parameters to your interactive job 
command. See :ref:`sge-interactive`.

* Do you need to run a parallel application on multiple processors or nodes?

If you have a task that uses distributed memory parallelism, like MPI, you
will need to read the section on ``mpie``.

* Do you want to visualise your results without downloading all the data from iceberg?

You can make use of iceberg's remote visualisation nodes, or run 
paraview or the Jupyter (IPython) notebook in server mode. See ``remote-vis``.


.. _running-applications:

Loading Applications or Libraries
---------------------------------

Once you have connected to iceberg you have to run `qsh` or `qrsh` to get a 
terminal on a woker node (see :ref:`ssh`).
Once you have a terminal on a worker node you can load software using the modules
system.
For information on modules and the applications availible see :ref:`software`.

