icu - International Components for Unicode
==========================================

.. sidebar:: icu

   :Latest Version: 58.2
   :URL: http://site.icu-project.org/

ICU is a mature, widely used set of C/C++ and Java libraries providing Unicode and Globalization support for software applications. ICU is widely portable and gives applications the same results on all platforms and between C/C++ and Java software.

ICU is released under a nonrestrictive open source license that is suitable for use with both commercial software and with other open source or free software.

Usage
-----
Version 58.2 of the icu library for C requires gcc version 4.9.2 (for the C++ standard library); To make the library **and this compiler** available, run the following: ::

        module load libs/icu/58.2/gcc-4.9.2

Installation Notes
------------------
This section is primarily for administrators of the system.

Version 58.2
^^^^^^^^^^^^

This Icu 58.2 build links against the GCC 4.9.2 C++ standard library and was installed as a dependency of `boost_sharc` (build using the same C++ standard library); Boost in turn was installed as a dependency of Caffe.

#. Download, configure, build, test and install using :download:`this script </sharc/software/install_scripts/icu/58.2/gcc-4.9.4/install.sh>`
#. Check the console output of the install process for ``All tests OK``
#. Install :download:`this modulefile </sharc/software/modulefiles/icu/58.2/gcc-4.9.4>` as ``/usr/local//modulefiles/icu/58.2/gcc-4.9.4``

