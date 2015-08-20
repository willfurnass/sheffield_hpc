Mathematica
===========

.. sidebar:: Wolfram Mathematica

   :Dependancies: None
   :URL: http://www.wolfram.com/mathematica/
   :Version: 10.2

Mathematica is a technical computing environment and programming language with strong symbolic and numerical abilities.

Single Core Interactive Usage
-----------------------------
After connecting to iceberg (see :ref:`ssh`),  start an interactive session with ``qsh``.

The latest version of Mathematica can be loaded with ::

        module load apps/binapps/mathematica

Alternatively, you can load a specific version of Mathematica using ::

        module load apps/binapps/mathematica/10.2

Mathematica can then be started with the ``mathematica`` command ::

        mathematica

Multicore Interactive Usage
---------------------------
Mathematica has extensive parallel functionality. To use it, you should request a parallel interactive session. For example, to request 4 cores ::

    qsh -pe openmp 4

Load and launch Mathematica ::

    module load apps/binapps/mathematica
    mathematica

In Mathematica, let's time how long it takes to the first 20 million primes on 1 CPU core ::

    AbsoluteTiming[primelist = Table[Prime[k], {k, 1, 20000000}];]

When I tried this, I got 78 seconds. Your results may vary greatly. Now, let's launch 4 ParallelKernels and redo the calculation in parallel ::

    LaunchKernels[4]
    AbsoluteTiming[primelist =
   ParallelTable[Prime[k], {k, 1, 20000000},
    Method -> "CoarsestGrained"];]

When I tried this, I got 29 seconds -- around 2.7 times faster. This illustrates a couple of points:-

* You should always request the same number of kernels as you requested in your ``qsh`` command (in this case, 4). If you request more, you will damage performance for yourself and other users of the system.
* N kernels doesn't always translate to N times faster.

Installation notes
------------------
These are primarily for administrators of the system ::

  mkdir -p /usr/local/packages6/apps/binapps/mathematica/10.2
  chmod +x ./Mathematica_10.2.0_LINUX.sh
  ./Mathematica_10.2.0_LINUX.sh

The installer is interactive. Here's the session output ::

  ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
                                                                                      Wolfram Mathematica 10.2 Installer
  ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

  Copyright (c) 1988-2015 Wolfram Research, Inc. All rights reserved.

  WARNING: Wolfram Mathematica is protected by copyright law and international treaties. Unauthorized reproduction or distribution may result in severe civil and criminal penalties and will be
  prosecuted to the maximum extent possible under law.

  Enter the installation directory, or press ENTER to select /usr/local/Wolfram/Mathematica/10.2:
  >

  Error: Cannot create directory /usr/local/Wolfram/Mathematica/10.2.

  You may need to be logged in as root to continue with this installation.

  Enter the installation directory, or press ENTER to select /usr/local/Wolfram/Mathematica/10.2:
  > /usr/local/packages6/apps/binapps/mathematica/10.2

  Now installing...

  [*********************************************************************************************************************************************************************************************************]

  Type the directory path in which the Wolfram Mathematica script(s) will be created, or press ENTER to select /usr/local/bin:
  > /usr/local/packages6/apps/binapps/mathematica/10.2/scripts

  Create directory (y/n)?
  > y


  WARNING: No Avahi Daemon was detected so some Kernel Discovery features will not be available. You can install Avahi Daemon using your distribution's package management system.

  For Red Hat based distributions, try running (as root):

  yum install avahi

  Installation complete.

Remove the ``playerpass`` file ::

  rm /usr/local/packages6/apps/binapps/mathematica/10.2/Configuration/Licensing/playerpass

Install the University network ``mathpass`` file at ``/usr/local/packages6/apps/binapps/mathematica/10.2/Configuration/Licensing``
