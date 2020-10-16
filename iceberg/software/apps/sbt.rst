sbt
===

.. sidebar:: sbt

   :Versions:  0.13.13
   :URL: http://www.scala-sbt.org/

sbt is a build tool for Scala, Java, and more.

Interactive Usage
-----------------
After connecting to Iceberg (see :ref:`ssh`),  start an interactive session with the :ref:`qrshx` command.

The latest version of sbt (currently 0.13.13) is made available with the command ::

        module load apps/binapps/sbt

Alternatively, you can load a specific version with ::

        module load apps/binapps/sbt/0.13.13

This command makes the sbt command available to your session.

Installation notes
------------------
These are primarily for administrators of the system.

sbt 0.13.13 was installed as follows ::

  tar -xvzf ./sbt-0.13.13.tgz
  cd sbt-launcher-packaging-0.13.13/
  mkdir -p /usr/local/packages6/apps/binapps/scala-sbt/0.13.13
  mv ./* /usr/local/packages6/apps/binapps/scala-sbt/0.13.13/

Module files
------------
* The module file is on the system at ``/usr/local/modulefiles/apps/binapps/sbt``
* On github: :download:`0.13.13 </iceberg/software/modulefiles/apps/binapps/sbt/0.13.13>`.
