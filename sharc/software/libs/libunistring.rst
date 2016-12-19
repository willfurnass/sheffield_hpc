libunistring
============

.. sidebar:: libunistring

   :Latest Version: 0.9.7
   :URL: http://www.gnu.org/software/libunistring/

Text files are nowadays usually encoded in Unicode, and may consist of very
different scripts - from Latin letters to Chinese Hanzi, with many kinds of
special characters: accents, right-to-left writing marks, hyphens, Roman
numbers, and much more. But the POSIX platform APIs for text do not contain
adequate functions for dealing with particular properties of many Unicode
characters. In fact, the POSIX APIs for text have several assumptions at their
base which don’t hold for Unicode text.

This library provides functions for manipulating Unicode strings and for
manipulating C strings according to the Unicode standard.

Usage
-----
To make the library available, run the following: ::

       module load libs/libunistring/0.9.7/gcc-4.9.4

This correctly populates the environment variables ``LD_LIBRARY_PATH``, ``LIBRARY_PATH`` and ``CPATH``.

Installation Notes
------------------
This section is primarily for administrators of the system.

Version 0.9.7
^^^^^^^^^^^^^

This build was installed as a dependency of `boost_sharc` (build using the same C++ standard library); Boost in turn was installed as a dependency of Caffe.

#. Download, configure, build, test and install using :download:`this script </sharc/software/install_scripts/libunistring/0.9.7/gcc-4.9.4/install.sh>`
#. Check the console output of the install process to check that no tests have errored/failed: ``TOTAL: 499 / PASS: 489 / SKIP: 10``
#. Install :download:`this modulefile </sharc/software/modulefiles/libunistring/0.9.7/gcc-4.9.4>` as ``/usr/local//modulefiles/libunistring/0.9.7/gcc-4.9.4``


