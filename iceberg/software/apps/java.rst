.. Warning:: 
    Iceberg reaches end-of-life on **30th November 2020.**
    If you are running jobs on Iceberg then you need to take urgent action to ensure that your jobs/scripts will run on ShARC or Bessemer. 
 
    If you have never used ShARC or Bessemer then now is the time to test your scripts.
    Not all software on Iceberg is available on ShARC/Bessemer. 

.. _Java-iceberg:

Java
====

.. sidebar:: Java

   :Latest Version: 1.8.0u112
   :URL: https://www.java.com/en/download/

Java is a programming language and computing platform first released by Sun Microsystems in 1995.

Interactive Usage
-----------------
After connecting to Iceberg (see :ref:`ssh`), start an interactive session with the `qrsh` or `qsh` command.

The latest version of Java (currently 1.8.0u112) is made available with the command: ::

        module load apps/java

Alternatively, you can load a specific version with one of ::

        module load apps/java/1.8.0u112
        module load apps/java/1.8u71
        module load apps/java/1.7.0u55
        module load apps/java/1.7
        module load apps/java/1.6

Check that you have the version you expect. First, the runtime: ::

        java -version

        java version "1.8.0_112"
        Java(TM) SE Runtime Environment (build 1.8.0_112-b15)
        Java HotSpot(TM) 64-Bit Server VM (build 25.112-b15, mixed mode)

Now, the compiler ::

        javac -version

        javac 1.8.0_112

Virtual Memory
--------------

.. note::

   The following is only relevant when revisiting older job submission scripts and documentation 
   as the job scheduler :ref:`now polices jobs using real memory, not virtual memory <real-vs-virt-mem>`.

By default, Java requests a lot of *virtual memory* on startup.
This is usually a given fraction of the *physical memory* on a node,
which can be quite a lot on Iceberg.
See :ref:`real-vs-virt-mem` for explanations of the difference between virtual and real/physical memory.

The default amount of virtual memory that Java uses at startup is controlled in two different ways on Iceberg:

For versions **prior to 1.8.0u112** we created a wrapper script to ``java`` that uses the ``-Xmx1G`` switch to force Java to restrict its *heap size* to a maximum of 1 GB of memory.  If this is amount is insufficient, you are free to allocate as much memory as you require by setting your own value for ``-Xmx``.

For example, consider the following submission script: ::

  #!/bin/bash
  #Request 9 gigabytes of real memory from the scheduler (mem)
  #$ -l rmem=9G

  # load the Java module
  module load apps/java/1.8u71

  #Run java program allocating 5 Gigabytes
  java -Xmx5G HelloWorld

Note that **the above is not possible** if an appliation (e.g. a startup script for another software package on the cluster) starts ``java`` itself (instead of you explicitly starting it)
and the above also will not help if an application uses an internally packaged version of Java (rather than one that can be *activated* using ``module load``).

For **versions 1.8.0u112 onwards** the maximum heap size is instead to be restricted using an *environment variable*.  The following is set when you run ``module load ...``: ::

        _JAVA_OPTIONS='-Xmx1G'

You can override this default value by running something like: ::

        export _JAVA_OPTIONS='-Xmx6G'

before starting your application that depends on Java.
``_JAVA_OPTIONS`` can be interpretted by Java programs you start and Java programs started by other programs,
as well as by Java Virtual Machines (JVMs) that you activate using ``module load`` and JVMs that are packaged within applications.

Installation notes
------------------

These are primarily for administrators of the system.

Version 1.8.0u112
^^^^^^^^^^^^^^^^^

#. Unpack the Java JDK tarball into ``/usr/local/packages6/apps/binapps/java/jdk1.8.0u112/``
#. Install :download:`this modulefile </iceberg/software/modulefiles/apps/binapps/java/1.8.0u112>` as ``/usr/local/modulefiles/apps/binapps/java/1.8.0u112``.

Note that the modulefile contains the following line: ::

        setenv _JAVA_OPTIONS -Xmx1G

Versions prior to 1.8.0u112
^^^^^^^^^^^^^^^^^^^^^^^^^^^

Unzip and copy the install directory to ``/usr/local/packages6/apps/binapps/java/jdk${VERSION}/``

To fix the virtual memory issue described above, we use a wrapper around the java install that sets Java's maximum heap size (``Xmx``) parameter to a reasonable value.

Create the file ``/usr/local/packages6/apps/binapps/java/jdk1.8.0_71/shef/java`` with contents: ::

        #!/bin/bash
        #
        # Java version 1.8 cannot be invoked without specifying the java virtual
        # machine size due to the limitations imposed by us via SGE on memory usage.
        # Therefore this script intercepts the java invocations and adds a
        # memory constraint parameter to java engine unless there was one already
        # specified on the command parameter.
        #
        #
          if test -z "`echo $* | grep -e -Xmx`"; then
        # user has not specified -Xmx memory requirement flag, so add it.
            /usr/local/packages6/apps/binapps/java/jdk1.8.0_71/bin/java -Xmx1G $*
        else
        # user specified the -Xmx flag, so don't add it.
            /usr/local/packages6/apps/binapps/java/jdk1.8.0_71/bin/java $*
        fi

The module file is at ``/usr/local/modulefiles/apps/java/1.8u71``. Its contents are ::

  #%Module10.2#####################################################################

  ## Module file logging
  source /usr/local/etc/module_logging.tcl
  ##


  proc ModulesHelp { } {
      global helpmsg
      puts stderr "\t$helpmsg\n"
  }


  set version 1.8

  set javahome /usr/local/packages6/apps/binapps/java/jdk1.8.0_71/

  if [ file isdirectory $javahome/bin ] {
      module-whatis "Sets JAVA to version $version"
      set helpmsg "Changes the default version of Java to Version $version"
      # bring in new version
      setenv JAVA_HOME $javahome
      prepend-path PATH $javahome/bin
      prepend-path PATH $javahome/shef
      prepend-path MANPATH $javahome/man
  } else {
      module-whatis "JAVA $version not installed"
      set helpmsg "JAVA $version not installed"
      if [ expr [ module-info mode load ] || [ module-info mode display ] ] {
  	# bring in new version
  	puts stderr "JAVA $version not installed on [uname nodename]"
      }
  }
