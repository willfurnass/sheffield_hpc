sbt
===

.. sidebar:: sbt

   :Versions:  0.13.13
   :URL: http://www.scala-sbt.org/

SBT is a build tool for Scala, Java, and more.

Interactive Usage
-----------------
After connecting to Iceberg (see :ref:`ssh`), :ref:`start an interactive session <sched_interactive>` then
load a specific version of SBT with ::

        module load apps/binapps/sbt/0.13.13

This command makes the ``sbt`` command available to your session.

Installation notes
------------------

These are primarily for administrators of the system.

SBT 0.13.13 was installed as follows ::

   tar -xvzf ./sbt-0.13.13.tgz
   cd sbt-launcher-packaging-0.13.13/
   mkdir -p /usr/local/packages6/apps/binapps/scala-sbt/0.13.13
   mv ./* /usr/local/packages6/apps/binapps/scala-sbt/0.13.13/

Module files
------------

:download:`This modulefile </iceberg/software/modulefiles/apps/binapps/sbt/0.13.13>`
is installed at ``/usr/local/modulefiles/apps/binapps/sbt``.
