.. Warning:: 
    Iceberg reaches end-of-life on **30th November 2020.**
    If you are running jobs on Iceberg then you need to take urgent action to ensure that your jobs/scripts will run on ShARC or Bessemer. 
 
    If you have never used ShARC or Bessemer then now is the time to test your scripts.
    Not all software on Iceberg is available on ShARC/Bessemer. 

htop
====

.. sidebar:: htop

   :Versions:  2.0
   :URL: http://hisham.hm/htop/

This is htop, an interactive process viewer for Unix systems.

Interactive Usage
-----------------
After connecting to Iceberg (see :ref:`ssh`),  start an interactive session with the `qrsh` or `qsh` command.

The latest version of htop (currently 2.0) is made available with the command ::

        module load apps/gcc/4.4.7/htop

Alternatively, you can load a specific version with ::

        module load apps/gcc/4.4.7/htop/2.0

This command makes the `htop` binary available to your session.

Installation notes
------------------
htop was installed using gcc 4.4.7 ::

    tar -xvzf ./htop-2.0.0.tar.gz
    cd htop-2.0.0
    mkdir -p /usr/local/packages6/apps/gcc/4.4.7/htop/2.0
    ./configure --prefix=/usr/local/packages6/apps/gcc/4.4.7/htop/2.0
    make
    make install

Testing
-------
No test suite was found.

Modulefile
----------
* The module file is on the system at ``/usr/local/modulefiles/apps/gcc/4.4.7/htop/2.0``
* The module file is :download:`on github </iceberg/software/modulefiles/apps/gcc/4.4.7/htop/2.0>`.
