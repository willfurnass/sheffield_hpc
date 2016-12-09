sbt
===

.. sidebar:: sbt

   :Versions:  0.13.13
   :URL: http://www.scala-sbt.org/

sbt is a build tool for Scala, Java, and more.

Interactive Usage
-----------------
After connecting to ShARC, start an interactive session with the :ref:`qrshx` command.

The latest version of sbt (currently 0.13.13) is made available with the command ::

        module load dev/sbt

Alternatively, you can load a specific version with ::

        module load dev/sbt/0.13.13

This command makes the `sbt` command available to your session.
You will also need to load a version of Java. For example ::

       module load apps/java/jdk1.8.0_102/binary

Installation notes
------------------
These are primarily for administrators of the system.

sbt 0.13.13 was installed as follows ::

  cp /usr/local/media/sbt/sbt-0.13.13.tgz .
  tar -xvzf ./sbt-0.13.13.tgz
  cd sbt-launcher-packaging-0.13.13/
  mkdir -p /usr/local/packages/dev/sbt/0.13.13
  mv ./* /usr/local/packages/dev/sbt/0.13.13

Installation Media
------------------
* The install tarball is on the system at ``/usr/local/media/sbt/sbt-0.13.13.tgz``

Module files
------------
* The module file is on the system at ``/usr/local/modulefiles/dev/sbt/0.13.13``
* On github: :download:`0.13.13 </sharc/software/modulefiles/dev/sbt/0.13.13>`.
