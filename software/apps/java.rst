Java
====

.. sidebar:: Java

   :Latest Version: 1.8u71
   :URL: https://www.java.com/en/download/

Java is a programming language and computing platform first released by Sun Microsystems in 1995.

Interactive Usage
-----------------
After connecting to Iceberg (see :ref:`ssh`), start an interactive session with the `qrsh` or `qsh` command. For Java, we suggest that you request a large amount of memory ::

    qrsh -l mem=32G -l rmem=32G

The latest version of Java (currently 1.8u71) is made available with the command ::

        module load apps/java

Alternatively, you can load a specific version with one of ::

       module load apps/java/1.8u71
       module load apps/java/1.7.0u55
       module load apps/java/1.7
       module load apps/java/1.6

Check that you have the version you expect. First, the runtime ::

    java -version

    java version "1.8.0_71"
    Java(TM) SE Runtime Environment (build 1.8.0_71-b15)
    Java HotSpot(TM) 64-Bit Server VM (build 25.71-b15, mixed mode)

Now, the compiler ::

    javac -version

    javac 1.8.0_71

Installation notes
------------------
These are primarily for administrators of the system.

The module file is at ``/usr/local/modulefiles/apps/java/1.8u71``. It's contents are ::

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
