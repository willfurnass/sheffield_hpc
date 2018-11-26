ffmpeg
======

.. sidebar:: ffmpeg

   :Version: 4.1
   :URL: https://www.ffmpeg.org/

FFmpeg is the leading multimedia framework, able to decode, encode, transcode, mux, demux, stream, filter and play pretty much anything that humans and machines have created.
It supports the most obscure ancient formats up to the cutting edge.
No matter if they were designed by some standards committee, the community or a corporation.
It is also highly portable: FFmpeg compiles, runs, and passes our testing infrastructure FATE across Linux, Mac OS X, Microsoft Windows, the BSDs, Solaris, etc. under a wide variety of build environments, machine architectures, and configurations.

Interactive Usage
-----------------
After connecting to iceberg (see :ref:`ssh`),  start an interactive session with the :code:`qsh` or :code:`qrsh` command.
The latest version of ffmpeg (currently 4.1) is made available with the command ::

        module load apps/ffmpeg

Alternatively, you can load a specific version with ::

        module load apps/ffmpeg/4.1/gcc-4.9.4

This command makes the ffmpeg binaries available to your session. It also loads version 4.9.4 of the gcc compiler environment since gcc 4.9.4 was used to compile ffmpeg 4.1

You can now run ffmpeg. For example, to confirm the version loaded ::

    ffmpeg -version

and to get help ::

    ffmpeg -h

Documentation
-------------
Once you have made ffmpeg available to the system using the `module` command above, you can read the man pages by typing ::

    man ffmpeg

Installation notes
------------------

LAMMPS was compiled using the
:download:`install_ffmpeg.sh </sharc/software/install_scripts/apps/ffmpeg/4.1/gcc-4.9.4/install_ffmpeg.sh>` script.

The module file is
:download:`/usr/local/modulefiles/apps/ffmpeg/4.1/gcc-4.9.4 </sharc/software/modulefiles/apps/ffmpeg/4.1/gcc-4.9.4>`.

Testing
-------
The test suite was executed (immediately post make install) ::

    make check

All tests passed.
