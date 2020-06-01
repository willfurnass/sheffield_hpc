.. _doxygen_sharc:

Doxygen
=======

Doxygen is a tool for building documentation for source code e.g. inter-related HTML pages for C++ source code.

Usage
-----

Doxygen can be loaded with: ::

    module load apps/doxygen/1.8.13/gcc-4.9.4

NB Doxygen has a run-time dependency on ``libstdc++`` so the above also needs to
(and does) load the :ref:`GCC compiler <gcc_sharc>` version 4.9.4.

Installation
------------

Version 1.8.13
^^^^^^^^^^^^^^

1. Install using :download:`this script </sharc/software/install_scripts/dev/doxygen/1.8.13/gcc-4.9.4/install.sh>`
2. Install :download:`this modulefile </sharc/software/modulefiles/dev/doxygen/1.8.13/gcc-4.9.4>` as ``/usr/local/modulefiles/dev/doxygen/1.8.13/gcc-4.9.4``
