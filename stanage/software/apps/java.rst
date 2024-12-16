.. _java_stanage:

Java (Temurin/OpenJDK)
======================

.. sidebar:: Java

   :Latest Version: 17.0.4
   :URL: https://adoptium.net/en-GB/temurin/releases/

Java is a programming language and computing platform first released by Sun Microsystems in 1995.

OpenJDK is an open-source implementation of Java.

Interactive Usage
-----------------

.. include:: /referenceinfo/imports/scheduler/SLURM/common_commands/srun_start_interactive_session_import_stanage.rst

You can then load a version of of Java using one of the following:

.. include:: /referenceinfo/imports/stanage/packages/java-ml-el7-icelake-znver-stanage.rst

NB ``Java/11.0.2`` is `OpenJDK <https://openjdk.org/>`__;
``Java/11.0.16`` and newer are `Eclipse Temurin <https://adoptium.net/en-GB/temurin/releases/>`__, which is based on OpenJDK.


Check that you have the version you expect. First, the runtime:

.. code-block:: console
   :emphasize-lines: 1
   
   $ java -version
   openjdk version "17.0.4" 2022-07-19
   OpenJDK Runtime Environment Temurin-17.0.4+8 (build 17.0.4+8)
   OpenJDK 64-Bit Server VM Temurin-17.0.4+8 (build 17.0.4+8, mixed mode, sharing)

Now, the compiler:

.. code-block:: console
   :emphasize-lines: 1

   $ javac -version
   javac 17.0.4

Installation notes
------------------
This section is primarily for administrators of the system. Java has been installed using the default Easybuild config files.

Build logs and test reports can be found in ``$EBROOTJAVA/easybuild`` with a given module loaded.
