.. _filestore:

Filestores
==========

Every HPC user has access to *up to* five different storage areas:

* :ref:`home_dir`: per-user backed-up storage
* :ref:`data_dir`: additional per-user backed-up storage (*not on Bessemer*)
* :ref:`fastdata_dir`: high-performance shared filesystem for temporary data - optimised for reading/writing large files from multiple nodes and threads simultaneously
* :ref:`shared_dir`: per-PI shared storage areas for project data - can be accessed from non-HPC machines too
* :ref:`scratch_dir`: per-node temporary storage - useful for reading/writing lots of small files within *one job*

The storage areas differ in terms of:

* the amount of space available;
* whether they are available from multiple nodes;
* whether they are shared between clusters;
* whether the underlying storage system is performant if reading/writing large files;
* whether the underlying storage system is performant if reading/writing small files;
* frequency of backup and the time that the data can be left there.

.. _home_dir:

Home directories
----------------
All users have a home directory on each system, some of which are shared between systems:

+----------+------------------------+------+----------------+-----------------------------------------------+-------------------------+
| System   | Path                   | Type | Quota per user | Shared between system login and worker nodes? | Shared between systems? |
+==========+========================+======+================+===============================================+=========================+
| Bessemer | ``/home/yourusername`` | NFS  | 100GB          | Yes                                           | No                      |
+----------+------------------------+------+----------------+-----------------------------------------------+-------------------------+
| ShARC    | ``/home/yourusername`` | NFS  | 10GB           | Yes                                           | ShARC + Iceberg         |
+----------+------------------------+------+----------------+-----------------------------------------------+-------------------------+
| Iceberg  | ``/home/yourusername`` | NFS  | 10GB           | Yes                                           | ShARC + Iceberg         |
+----------+------------------------+------+----------------+-----------------------------------------------+-------------------------+

See also: :ref:`quota_check` and * :ref:`exceed_quota`.

Backups
^^^^^^^

+---------------------------+--------------------+---------------------------------------+
| Frequency of snapshotting | Snapshots retained | Mirrored onto separate storage system |
+===========================+====================+=======================================+
| Every 4 hours             | 10 most recent     | Yes                                   |
+---------------------------+--------------------+---------------------------------------+
| Every night               | Last 28 days       | Yes                                   |
+---------------------------+--------------------+---------------------------------------+

See also: :ref:`recovering_snapshots`.

.. _data_dir:

*Data* directories
------------------

Every user on Iceberg and ShARC (**not Bessemer**) has access to a larger *data* storage area:

+----------+------------------------+------+----------------+-----------------------------------------------+-------------------------+
| System   | Path                   | Type | Quota per user | Shared between system login and worker nodes? | Shared between systems? |
+==========+========================+======+================+===============================================+=========================+
| Bessemer | N/A                    | NFS  | N/A            | N/A                                           | N/A                     |
+----------+------------------------+------+----------------+-----------------------------------------------+-------------------------+
| ShARC    | ``/data/yourusername`` | NFS  | 100GB          | Yes                                           | ShARC + Iceberg         |
+----------+------------------------+------+----------------+-----------------------------------------------+-------------------------+
| Iceberg  | ``/data/yourusername`` | NFS  | 100GB          | Yes                                           | ShARC + Iceberg         |
+----------+------------------------+------+----------------+-----------------------------------------------+-------------------------+

See also: :ref:`quota_check` and * :ref:`exceed_quota`.

Backups
^^^^^^^

+---------------------------+--------------------+---------------------------------------+
| Frequency of snapshotting | Snapshots retained | Mirrored onto separate storage system |
+===========================+====================+=======================================+
| Every 4 hours             | 10 most recent     | No                                    |
+---------------------------+--------------------+---------------------------------------+
| Every night               | Last 7 days        | No                                    |
+---------------------------+--------------------+---------------------------------------+

See also: :ref:`recovering_snapshots`.

Automounting
^^^^^^^^^^^^^

*Data* directories are **made available to you (mounted) on demand**: 
if you list the contents of just ``/data`` after first logging on then ``/data/yourusername`` subdirectories might not be shown.
However, if you list the contents of ``/data/yourusername`` itself or change into that directory
then its contents will appear.  

Later on if you list the contents of ``/data`` again 
you may find that ``/data/yourusername`` has disappeared again, as 
it is automatically *unmounted* following a period of inactivity.  

.. _fastdata_dir:

*Fastdata* areas
----------------

**Fastdata** areas are **optimised for large file operations**.  
These areas are `Lustre <https://en.wikipedia.org/wiki/Lustre_(file_system)>`__ filesystems. 

They are are **faster** than :ref:`home_dir`, :ref:`data_dir` and :ref:`shared_dir` when dealing with larger files but 
are **not performant when reading/writing lots of small files** 
(:ref:`scratch_dir` are ideal for reading/writing lots of small temporary files within jobs).
An example of how slow it can be for large numbers of small files is detailed `here <http://www.walkingrandomly.com/?p=6167>`__.

+----------+---------------------+--------+----------------+---------------------+--------------------------------------------------------+---------------------------+
| System   | Path                | Type   | Quota per user | Filesystem capacity | Shared between systems?                                | Network bandwith per link |
+==========+=====================+========+================+=====================+========================================================+===========================+
| Bessemer | ``/fastdata``       | Lustre | None           | 460 TB              | No                                                     | 25Gb/s Ethernet           |
+----------+---------------------+--------+----------------+---------------------+--------------------------------------------------------+---------------------------+
| ShARC    | ``/fastdata``       | Lustre | None           | 669 TB              | ShARC's fastdata filesystem is accessible from Iceberg | 100Gb/s (*Omni-Path*)     |
+----------+---------------------+--------+----------------+---------------------+--------------------------------------------------------+---------------------------+
| Iceberg  | ``/fastdata-sharc`` | Lustre | None           | 669 TB              | ShARC's fastdata filesystem is accessible from Iceberg | 1Gb/s Ethernet            |
+----------+---------------------+--------+----------------+---------------------+--------------------------------------------------------+---------------------------+

Backups
^^^^^^^

Fastdata areas are **not backed up**.

Managing your files in fastdata areas
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

In order to avoid interference from other users' files 
it is **important** that you store your files in a directory created and named the same as your username. e.g. ::

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

Automatic file deletion
^^^^^^^^^^^^^^^^^^^^^^^

.. warning::

    **There are no quota controls in fastdata areas** but 
    **older files** are **automatically deleted**: 
    a report of files older than 60 days is regularly generated, 
    the owners of these files are then notified by email then 
    a week after the email(s) are sent the identified files are deleted. 

    We reserve the right to change this policy without warning in order to ensure efficient running of the service.

    It is important to therefore not use *fastdata* areas for long-term storage and 
    **copy important data** from these areas to **backed-up areas** (:ref:`home_dir`, :ref:`data_dir` or :ref:`shared_dir`).

You can use the ``lfs``  command to find out which files in a *fastdata* directory are older than a certain number of days and hence approaching the time of deletion. 
For example, to find files 50 or more days old ::

    lfs find -ctime +50 /fastdata/yourusername

File locking
^^^^^^^^^^^^

POSIX file locking is not enabled on these Lustre filesystems, 
which can cause issues for certain applications that require/expect it
(e.g. programs that create/use SQLite databases).

.. _shared_dir:

*Shared* (project) directories
------------------------------

Each PI at the University is entitled to request a `free 10 TB storage area for sharing data with their group and collaborators <https://sheffield.ac.uk/it-services/research-storage/using-research-storage>`__.
The capacity per area can be extended and additional shared areas can be purchased (both at a cost).

After one of these project storage areas has been requested/purchased it can be accessed in two ways:

* as a Windows-style (SMB) file share on machines other than ShARC/Iceberg using ``\\uosfstore.shef.ac.uk\shared\``;
* as a subdirectory of ``/shared`` on ShARC/Iceberg (you need to **explicitly request HPC access when you order storage from IT Services**).

Backups
^^^^^^^

+---------------------------+--------------------+---------------------------------------+
| Frequency of snapshotting | Snapshots retained | Mirrored onto separate storage system |
+===========================+====================+=======================================+
| Every 4 hours             | 10 most recent     | Yes                                   |
+---------------------------+--------------------+---------------------------------------+
| Every night               | Last 7 days        | Yes                                   |
+---------------------------+--------------------+---------------------------------------+

See also: :ref:`recovering_snapshots`.
  
Automounting
^^^^^^^^^^^^

Similar to :ref:`data_dir`, subdirectories beneath ``/shared`` are **mounted on demand** on the HPC systems: 
they may not be visible if you simply list the contents of the ``/shared`` directory but 
will be accessible if you ``cd`` (change directory) into a subdirectory e.g. ``cd /shared/my_group_file_share1``.

Specifics for Bessemer
^^^^^^^^^^^^^^^^^^^^^^

If you need to access a ``/shared`` area on Bessemer please contact `helpdesk@sheffield.ac.uk <helpdesk@sheffield.ac.uk>`__ to arrange this.


.. warning::

        * If you access a ``/shared`` directory stored in Sheffield from Bessemer then you may experience slower performance, espeicially for small files.
        * Network traffic between Bessemer and Sheffield Research Filestore is not encrypted when travelling between Sheffield and Leeds over JANET
        * ``/shared`` areas can be created on Bessemer's filestore system if you need faster access from Bessemer


Permissions behaviour
^^^^^^^^^^^^^^^^^^^^^

You may encounter strange permissions issues when running programs on HPC against the ``/shared`` areas 
e.g. ``chmod +x /shared/mygroup1/myprogram.sh`` fails.
Here we try to explain why.

Behind the scenes, the file server that provides this shared storage manages permissions using 
Windows-style `ACLs <https://en.wikipedia.org/wiki/Access_control_list>`_ 
(which can be set by area owners via the `Research Storage management web interface <https://sheffield.ac.uk/storage>`__.
However, the filesystem is mounted on a Linux cluster using NFSv4 so the file server therefore requires 
a means for mapping Windows-style permissions to Linux ones.  
An effect of this is that the Linux `mode bits <https://en.wikipedia.org/wiki/Modes_(Unix)>`_ for files/directories under ``/shared`` on the HPC systems
are not always to be believed: 
the output of ``ls -l somefile.sh`` may indicate that a file is readable/writable/executable when 
the ACLs are what really determine access permissions.  
Most applications have robust ways of checking for properties such as executability but 
some applications can cause problems when accessing files/directories on ``/shared`` by naively checking permissions just using Linux mode bits:

* `which <http://linux.die.net/man/1/which>`_: 
  a directory under ``/shared`` may be on your path and 
  you may be able to run a contained executable without prefixing it with a absolute/relative directory 
  but ``which`` may fail to find that executable.
* Perl: scripts that check for executability of files on ``/shared`` using ``-x`` may fail 
  unless Perl is explicitly told to test for file permissions in a more thorough way 
  (see the mention of ``use filetest 'access'`` `here <http://perldoc.perl.org/functions/-X.html>`_).
* git: may complain that permissions have changed if 
  a repository is simply moved to ``/shared/someplace`` from elsewhere on Bessemer/ShARC/Iceberg. 
  As a workaround you can tell git to not to track Linux permissions for a single repository using 
  ``git config core.filemode false`` or 
  for all repositories using ``git config --global core.filemode false``.

**Changing how attempts to change permissions are handled**: each ``/shared`` area can be configured so that

#. Attempts to change file/directory mode bits fail (e.g. ``chmod +x /shared/mygroup1/myprogram.sh`` fails) (**default configuration per area**) **or**
#. Attempts to change file/directory mode bits appear to succeed (e.g. ``chmod +x /shared/mygroup1/myprogram.sh`` does not fail but also does not actually change any permissions on the underlying file server) (**alternative configuration per area**)

Contact the Helpdesk if you would like to switch to using the second way of handling permissions for a particular ``/shared/`` area.

Further information
^^^^^^^^^^^^^^^^^^^

The documentation for the ``/shared`` storage service includes information on:

* `how access/permissions are managed <https://www.sheffield.ac.uk/it-services/research-storage/access-rights>`__
* `how to create folders with associated permissions <https://www.sheffield.ac.uk/it-services/research-storage/create-folders>`__ 
  within ``/shared`` storage areas

.. _scratch_dir:

*Scratch* directories
---------------------

For **jobs that need to read/write lots of small files** the most performant storage will be 
the temporary storage on each node (under the ``/scratch`` directory).

This is because with :ref:`home_dir`, :ref:`data_dir`, :ref:`fastdata_dir`, :ref:`shared_dir`
each time a file is accessed the filesystem needs to request ownership/permissions information from another server
and for small files these overheads are proportionally high. 
However, for ``/scratch`` such ownership/permissions metadata is available on the local machine, 
so it is faster when dealing with small files.

The most obvious disadvantage to the ``/scratch`` node-local storage is that 
a given directory cannot reliabily be accessed between jobs as
you cannot guarantee that your next job will run on the same node.
Any data of value must therefore be **copied off** ``/scratch`` 
(e.g. to :ref:`home_dir` or :ref:`data_dir`)
**before the end of your job**.

**Where to store data beneath** ``/scratch``: 
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


.. _quota_check:

How to check your quota usage
-----------------------------

To find out your storage quota usage for your :ref:`home directory <home_dir>`, :ref:`data directory <data_dir>` (if not on Bessemer) and particular :ref:`shared_dir`: ::

    df -h somedirectoryname

For example:

+--------------------------------------------------+------------------------------+
| Storage area                                     | Command to check quota       |
+==================================================+==============================+
| :ref:`Home directory <home_dir>`                 | ``df -h /home/$USER``        |
+--------------------------------------------------+------------------------------+
| :ref:`Data directory <data_dir>`                 | ``df -h /data/$USER``        |
+--------------------------------------------------+------------------------------+
| A :ref:`Shared (project) directory <shared_dir>` | ``df -h /shared/myproject1`` |
+--------------------------------------------------+------------------------------+

.. _exceed_quota:

If you exceed your filesystem quota
-----------------------------------

If you reach your quota for your :ref:`home directory <home_dir>` then
many common programs/commands may cease to work as expected or at all and
you may not be able to log in.

In addition, jobs may fail if you exceed your quota
for your :ref:`data directory <data_dir>` or a :ref:`Shared (project) directory <shared_dir>`.

In order to avoid this situation it is strongly recommended that you:

* :ref:`Check your quota usage <quota_check>` regularly.
* Copy files that do not need to be backed up to a :ref:`Fastdata area <fastdata_dir>`
  or remove them from Bessemer/ShARC/Iceberg completely.

.. _recovering_snapshots:

Recovering files from backups
-----------------------------

:ref:`home_dir`, :ref:`data_dir` and :ref:`shared_dir` are regularly backed up.
See above for details of the backup schedules per area.
These backup processes create a series of storage area *snapshots*,
a subset of which can be accessed by HPC users from the HPC systems themselves
by *explicitly* browsing to hidden directories e.g.

+--------------------------------------------------+----------------------------------+
| Storage area                                     | Parent directory of snapshots    |
+==================================================+==================================+
| :ref:`Home directory <home_dir>`                 | ``/home/$USER/.snapshot``        |
+--------------------------------------------------+----------------------------------+
| :ref:`Data directory <data_dir>`                 | ``/data/$USER/.snapshot``        |
+--------------------------------------------------+----------------------------------+
| A :ref:`Shared (project) directory <shared_dir>` | ``/shared/myproject1/.snapshot`` |
+--------------------------------------------------+----------------------------------+

From within per-snapshot directories you can access (read-only) copies of files/directories.
This allows you to attempt recover any files you might have accidentally modified or deleted recently.

Note that ``.snapshot`` directories are not visible when listing all hidden items within their parent directories
(e.g. using ``ls -a /home/$USER``): 
you need to explicitly ``cd`` into ``.snapshot`` directories to see/access them.
