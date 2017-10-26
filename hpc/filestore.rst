.. _filestore:

Filestores
==========

Every user on ShARC/Iceberg has access to five different types of filestore:

* ``/home/yourusername``
* ``/data/yourusername``
* ``/fastdata``
* ``/shared/volumename``
* ``/scratch``

They differ in terms of:

* the amount of space available;
* whether they are available from multiple nodes;
* whether they are shared between clusters;
* whether the underlying storage system is performant if reading/writing large files;
* whether the underlying storage system is performant if reading/writing small files;
* frequency of backup and the time that the data can be left there.

Here are the current details of filestore available to each user.

``/home`` directory
-------------------
All users have a home directory in the location ``/home/yourusername``. 

The filestore quota is **10 GB** per user.

This area is shared between ShARC and Iceberg
and is accessible to all worker and login nodes.

**Backup policy:** ``/home`` has backup snapshots taken every 4 hours and 
we keep the 10 most recent. 
``/home`` also has daily snapshots taken each night, 
and we keep 28 days worth, 
mirrored onto a separate storage system.

The filesystem is NFS.

.. _data_dir:

``/data`` directory
-------------------
Every user has access to a much larger data-storage area provided at the location ``/data/yourusername``.

The quota for this area is **100 GB** per user.

This area is shared between ShARC and Iceberg 
and is accessible to all worker and login nodes.

**Backup policy:** ``/data`` has snapshots taken every 4 hours and we keep the 10 most recent. 
``/data`` also has daily snapshots taken each night, 
and we keep 7 days worth, 
but this is not mirrored.

The filesystem is NFS.

**Note**: the directory ``/data/yourusername`` is **made available to you (mounted) on demand**: 
if you list the contents of ``/data`` after first logging on then this subdirectory might not be shown.
However, if you list the contents of ``/data/yourusername`` itself or change into that directory
then its contents will appear.  
Later on if you list the contents of ``/data`` again 
you may find that ``/data/yourusername`` has disappeared again, as 
it is automatically *unmounted* following a period of inactivity.  

``/fastdata`` directory
-----------------------

All users also have access to a large fast-access data storage area under ``/fastdata``.  
This is **not** shared between ShARC and Iceberg.

In order to avoid interference from other users' files 
it is **vitally important** that you store your files in a directory created and named the same as your username. e.g. ::

    mkdir /fastdata/yourusername

By default the directory you create will have world-read access.  
If you want to restrict read access to just your account then run ::

    chmod 700 /fastdata/yourusername

after creating the directory. 
A more sophisticated sharing scheme would have private and public directories ::

    mkdir /fastdata/yourusername
    mkdir /fastdata/yourusername/public
    mkdir /fastdata/yourusername/private

    chmod 755 /fastdata/yourusername
    chmod 755 /fastdata/yourusername/public
    chmod 700 /fastdata/yourusername/private

The fastdata area provides **669 TeraBytes** of storage on SHARC and **260 Terabytes** of storage on Iceberg.
It takes advantage of the internal high-speed network infrastructure (OmniPath interconnects for ShARC; Infiniband for Iceberg) for fast access to data.

``/fastdata`` is **optimised for large file operations**.
it is faster than ``/home/``, ``/data`` and ``/shared`` (see below) when dealing with larger files but 
does **not handle lots of small files very well**:  it is less efficient than ``/scratch`` (see below) when dealing with smaller files.
An example of how slow it can be for large numbers of small files is detailed at http://www.walkingrandomly.com/?p=6167

Iceberg users: although ``/fastdata`` is available on all the worker nodes, 
you need to run your jobs on the newer nodes with Intel CPUs for best performance (by specifying ``-l arch=intel*`` in your job submission script).

.. warning::

    **There are no quota controls on the** ``/fastdata`` **area** but 
    **older files** on /fastdata are **automatically deleted**: 
    a report of files older than 60 days is regularly generated, 
    the owners of these files are then notified by email then 
    a week after the email(s) are sent the identified files are deleted. 

    We reserve the right to change this policy without warning in order to ensure efficient running of the service.

    It is important to therefore not use /fastdata for long-term storage and 
    **copy important data** on ``/fastdata`` to **backed-up areas** such as ``/home``, ``/data`` or ``/shared``.

You can use the ``lfs``  command to find out which files under ``/fastdata`` are older than a certain number of days and hence approaching the time of deletion. 
For example, to find files 50 or more days old ::

    lfs find -ctime +50 /fastdata/yourusername

**Backup policy:** ``/fastdata`` is **not backed up**.

``/fastdata`` uses the `Lustre <https://en.wikipedia.org/wiki/Lustre_(file_system)>`__ filesystem. 
This does not support POSIX locking which can cause issues for some applications 
(e.g. programs that create/use SQLite databases).

``/shared`` directories
-----------------------

CiCS now provide `10 terabytes of shared storage for free per research group <shef.ac.uk/cics/research-storage/using-research-storage>`__.
After the storage has been requested/purchased by a group's PI and then provisioned by CiCS it can be accessed by name

* as a Windows-style (SMB) file share on machines other than ShARC/Iceberg using ``\\uosfstore.shef.ac.uk\shared\``;
* as a subdirectory of ``/shared`` on ShARC/Iceberg (you need to **explicitly request HPC access when you order storage from CiCS**).
  
Note that this subdirectory will be **mounted on demand** on ShARC/Iceberg: 
it will not be visible if you simply list the contents of the ``/shared`` directory but 
will be accessible if you ``cd`` (change directory) into it e.g. ``cd /shared/my_group_file_share1``

**Regarding permissions**: 
behind the scenes, the file server that provides this shared storage manages permissions using 
Windows-style `ACLs <https://en.wikipedia.org/wiki/Access_control_list>`_ 
(which can be set by area owners via the `Research Storage management web interface <sheffield.ac.uk/storage>`__.
However, the filesystem is mounted on a Linux cluster using NFSv4 so the file server therefore requires 
a means for mapping Windows-style permissions to Linux ones.  
An effect of this is that the Linux `mode bits <https://en.wikipedia.org/wiki/Modes_(Unix)>`_ as seen on ShARC/Iceberg 
are not always to be believed for files under ``/shared``: 
the output of ``ls -l somefile.sh`` may indicate that a file is readable/writable/executable when 
the ACLs are what really determine access permissions.  
Most applications have robust ways of checking for properties such as executability but 
some applications can cause problems when accessing files/directories on ``/shared`` by naievely checking permissions just using Linux mode bits:

* `which <http://linux.die.net/man/1/which>`_: 
  a directory under ``/shared`` may be on your path and 
  you may be able to run a contained executable without prefixing it with a absolute/relative directory 
  but ``which`` may fail to find that executable.
* Perl: scripts that check for executability of files on ``/shared`` using ``-x`` may fail 
  unless Perl is explicitly told to test for file permissions in a more thorough way 
  (see the mention of ``use filetest 'access'`` `here <http://perldoc.perl.org/functions/-X.html>`_).
* git: may complain that permissions have changed if 
  a repository is simply moved to ``/shared/someplace`` from elsewhere on ShARC/Iceberg.  
  As a workaround you can tell git to not to track Linux permissions for a single repository using 
  ``git config core.filemode false`` or 
  for all repositories using ``git config --global core.filemode false``.

The documentation for the ``/shared`` storage serivce includes information on:

* `how access/permissions are managed <https://www.sheffield.ac.uk/cics/research-storage/access-rights>`__
* `how to create folders with associated permissions <https://www.sheffield.ac.uk/cics/research-storage/create-folders>`__ 
  within ``/shared`` storage areas

``/scratch``: for reading/writing small files
---------------------------------------------

For **jobs that need to read/write lots of small files** the most performant storage will be 
the temporary storage on each node (under the ``/scratch`` directory).

This is because with ``/home``, ``/data``, ``/fastdata`` and ``/shared`` 
each time a file is accessed the filesystem needs to request ownership/permissions information from another server
and for small files these overheads are proportionally high. 
However, for ``/scratch`` such ownership/permissions metadata is available on the local machine, 
so it is faster when dealing with small files.

The most obvious disadvantage to the ``/scratch`` node-local storage is that 
a given directory cannot relabily be accessed between jobs as
you cannot guarantee that your next job will run on the same node.
Any data of value must therefore be **copied off** ``/scratch`` 
(e.g. to ``/home`` or ``/data``)
**before the end of your job**.

**Where to store data within ``/scratch``**: 
The scheduler automatically creates a per-job directory for you under ``/scratch``.
If you started your job using ``qrshx``, ``qsh`` or ``qsub`` then 
the name of this directory is stored in the ``$TMPDIR`` environment variable e.g. ::

    [te1st@sharc-login1 ~]$ qrshx
    [te1st@sharc-node003 ~]$ cd $TMPDIR
    [te1st@sharc-node003 667443.1.all.q]$ pwd
    /scratch/667443.1.all.q

The scheduler will then clean up (delete) ``$TMPDIR`` at the end of your job, 
ensuring that the space can be used by other users.

If using ``qrsh`` to start your job then the environment variable will unfortunately be undefined
so you will need to manually create a directory under ``/scratch`` (named using your username)
and this will not be cleaned up when the job ends.

Anything under the ``/scratch`` may be deleted periodically when the worker-node is idle. 
``/scratch`` is **not backed up**.  There are no quotas for ``/scratch`` storage.

``/scratch`` uses the ext4 filesystem.

Determining your current filestore allocation
---------------------------------------------

To find out your current storage quota usage for ``/home`` and ``/data``: ::

    quota

If you exceed your file storage allocation
------------------------------------------

As soon as the quota is exceeded your account becomes frozen. 
In order to avoid this situation it is strongly recommended that you:

* Use the ``quota`` command to check your usage regularly.
* Copy files that do not need to be backed up to the  ``/fastdata/username`` area, 
  or remove them from ShARC/Iceberg completely.

Recovering snapshots 
--------------------

We take regular back-ups of your ``/home`` and ``/data`` directories and it is possible to directly access a limited subset of them.

There are 7 days worth of snapshots available in your ``/home`` and ``/data`` directories in 
a hidden directory called ``.snapshot``. 
You need to explicitly ``cd`` into this directory to get at the files::

    cd /home/YOURUSERNAME/.snapshot

The files are read-only. 
This allows you to attempt recover any files you might have accidentally deleted recently.

This does not apply for ``/fastdata`` for which we take no back-ups.
