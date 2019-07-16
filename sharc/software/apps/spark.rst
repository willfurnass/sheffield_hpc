.. _sparc_sharc:

spark
=====

.. sidebar:: Spark

   :Version: 2.3.0, 2.1
   :URL: http://spark.apache.org/

Apache Spark is a fast and general engine for large-scale data processing.

Interactive Usage
-----------------
After connecting to ShARC (see :ref:`ssh`),  start an interactive session with the :ref:`qrsh` or :ref:`qrshx` command.

Before using Spark, you will need to load a :ref:`version of Java <Java-sharc>`. For example: ::

    module load apps/java/jdk1.8.0_102/binary

To make a version of Spark available, use one the following commands: ::

    module load apps/spark/2.3.0/jdk-1.8.0_102
    module load apps/spark/2.1.0/gcc-4.8.5

You can now start a Spark shell session with ::

    spark-shell

SparkR
------
To use SparkR, you will additionally need to load a version of R e.g.: ::

    module load apps/java/jdk1.8.0_102/binary
    module load apps/spark/2.3.0/jdk-1.8.0_102
    module load apps/R/3.3.2/gcc-4.8.5

Now you can start a SparkR session by running: ::

    sparkR

Setting the number of cores
---------------------------
The installation of Spark on ShARC is limited to jobs that make use of one node.
As such, the maximum number of CPU cores you can request for a Spark job is (typically) 16.

First, you must request cores from the scheduler.
That is, you add the following line to your submission script to request 4 cores ::

  #$ -pe smp 4

You must also tell Spark to only use 4 cores by setting the ``MASTER`` environment variable ::

  export MASTER=local[4]

A full example using Python is given `here <https://github.com/mikecroucher/HPC_Examples/tree/master/languages/Python/pyspark_pi>`__.

Installation notes
------------------
These notes are primarily for administrators of the system.

Spark 2.3.0
^^^^^^^^^^^

* Install script: :download:`install.sh </sharc/software/install_scripts/apps/spark/2.3.0/jdk-1.8.0_102/install.sh>`
* Module file :download:`apps/spark/2.3.0/jdk-1.8.0_102 </sharc/software/modulefiles/apps/spark/2.3.0/jdk-1.8.0_102>`,
  which 

  * sets ``SPARK_HOME``
  * prepends the Spark ``bin`` directory to the ``PATH``
  * sets ``MASTER`` to ``local\[1\]`` (i.e. Spark will default to using 1 core)

Spark 2.1
^^^^^^^^^

.. code-block:: sh

   qrsh -l rmem=10G

   module load apps/java/jdk1.8.0_102/binary
   tar -xvzf ./spark-2.1.0.tgz
   cd spark-2.1.0
   ./build/mvn -DskipTests clean package

   mkdir -p /usr/local/packages/apps/spark/2.1
   cd ..
   mv spark-2.1.0 /usr/local/packages/apps/spark/2.1

The default install of Spark is incredibly verbose. Even a 'Hello World' program results in many lines of ``[INFO]``.
To make it a little quieter, the default log4j level has been reduced from ``INFO`` to ``WARN``: ::

    cd /usr/local/packages/apps/spark/2.1/spark-2.1.0/conf/
    cp log4j.properties.template log4j.properties
    
The file ``log4j.properties`` was then edited so that the line beginning ``log4j.rootCategory`` reads: ::
 
     log4j.rootCategory=WARN, console
     
Module file :download:`apps/spark/2.1/gcc-4.8.5 </sharc/software/modulefiles/apps/spark/2.1/gcc-4.8.5>`,
which 

* sets ``SPARK_HOME``
* prepends the Spark ``bin`` directory to the ``PATH``
* sets ``MASTER`` to ``local\[1\]`` (i.e. Spark will default to using 1 core)

