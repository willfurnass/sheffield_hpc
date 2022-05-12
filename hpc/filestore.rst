.. _filestore:

Filestores
==========

Every HPC user has access to *up to* six different storage areas:

* :ref:`home_dir`: per-user :term:`backed-up <Mirrored backups>`, :term:`snapshotted <Snapshotted storage>` storage
* :ref:`data_dir`: additional per-user snapshotted storage (*not on Bessemer*)
* :ref:`fastdata_dir`: high-performance shared filesystem for temporary data - optimised for reading/writing large files from multiple nodes and threads simultaneously
* :ref:`shared_dir`: per-PI shared storage areas (snapshotted and backed-up) for project data - can be accessed from non-HPC machines too
* :ref:`scratch_dir`: per-node temporary storage - useful for reading/writing lots of small files within *one job*
* :ref:`community_dir`: cluster-wide storage areas to allow users to share software.

The storage areas differ in terms of:

* the amount of space available;
* whether they are available from multiple nodes;
* whether they are shared between clusters;
* whether the underlying storage system is performant if reading/writing large files;
* whether the underlying storage system is performant if reading/writing small files;
* frequency of :term:`storage snapshotting <Snapshotted storage>`, 
  whether storage is :term:`mirrored <Mirrored backups>` 
  and the maximum duration data can be retained for;
* whether they handle permissions like a typical Linux filesystem.

-----

Choosing the correct filestore
------------------------------

To make a quick assessment of what storage area is likely to best fulfil your needs, please take a look at the provided decision tree below:

.. warning::

  This decision tree only provides a quick assessment, please check the full details of each filestore before committing to using them for your work.

.. raw:: html

  <div class="mxgraph" style="max-width:100%;border:1px solid transparent;" data-mxgraph="{&quot;highlight&quot;:&quot;#0000ff&quot;,&quot;nav&quot;:true,&quot;zoom&quot;:0.75,&quot;resize&quot;:true,&quot;toolbar&quot;:&quot;zoom layers tags lightbox&quot;,&quot;edit&quot;:&quot;_blank&quot;,&quot;xml&quot;:&quot;&lt;mxfile host=\&quot;app.diagrams.net\&quot; modified=\&quot;2022-05-12T13:11:51.554Z\&quot; agent=\&quot;5.0 (X11)\&quot; etag=\&quot;QhXgf6g9HzIt2syvAibr\&quot; version=\&quot;18.0.2\&quot; type=\&quot;google\&quot;&gt;&lt;diagram name=\&quot;Page-1\&quot; id=\&quot;9c096ad6-e400-ecc8-3e38-643d2caac077\&quot;&gt;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&lt;/diagram&gt;&lt;/mxfile&gt;&quot;}"></div>
  <script type="text/javascript" src="https://viewer.diagrams.net/js/viewer-static.min.js"></script>
  <!--This currently uses the remote JS hosted at diagrams.net but I saved a copy in the repo at _static/js/viewer-static.min.js in case service is ever withdrawn. -->


-----

.. _home_dir:

Home directories
----------------
All users have a home directory on each system, some of which are shared between systems:

+----------+------------------------+------+----------------+-----------------------------------------------+-------------------------+
| System   | Path                   | Type | Quota per user | Shared between system login and worker nodes? | Shared between systems? |
+==========+========================+======+================+===============================================+=========================+
| Bessemer | ``/home/yourusername`` | NFS  | 100GB          | Yes                                           | No                      |
+----------+------------------------+------+----------------+-----------------------------------------------+-------------------------+
| ShARC    | ``/home/yourusername`` | NFS  | 10GB           | Yes                                           | No                      |
+----------+------------------------+------+----------------+-----------------------------------------------+-------------------------+

See also: :ref:`quota_check` and * :ref:`exceed_quota`.

Snapshotting and mirrored backups
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

+---------------------------+--------------------+----------------------------------------+
| Frequency of snapshotting | Snapshots retained | Backed up onto separate storage system |
+===========================+====================+========================================+
| Every 4 hours             | 10 most recent     | Yes                                    |
+---------------------------+--------------------+----------------------------------------+
| Every night               | Last 28 days       | Yes                                    |
+---------------------------+--------------------+----------------------------------------+

See also: :ref:`recovering_snapshots`.

.. _data_dir:

*Data* directories
------------------

Every user on ShARC (**not Bessemer**) has access to a larger *data* storage area:

+----------+------------------------+------+----------------+-----------------------------------------------+-------------------------+
| System   | Path                   | Type | Quota per user | Shared between system login and worker nodes? | Shared between systems? |
+==========+========================+======+================+===============================================+=========================+
| Bessemer | N/A                    | NFS  | N/A            | N/A                                           | N/A                     |
+----------+------------------------+------+----------------+-----------------------------------------------+-------------------------+
| ShARC    | ``/data/yourusername`` | NFS  | 100GB          | Yes                                           | No                      |
+----------+------------------------+------+----------------+-----------------------------------------------+-------------------------+

See also: :ref:`quota_check` and * :ref:`exceed_quota`.

Snapshotting and mirrored backups
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

+---------------------------+--------------------+----------------------------------------+
| Frequency of snapshotting | Snapshots retained | Backed up onto separate storage system |
+===========================+====================+========================================+
| Every 4 hours             | 10 most recent     | No                                     |
+---------------------------+--------------------+----------------------------------------+
| Every night               | Last 7 days        | No                                     |
+---------------------------+--------------------+----------------------------------------+

See also: :ref:`recovering_snapshots`.

Automounting
^^^^^^^^^^^^^

*Data* directories are **made available to you (mounted) on demand**: 
if you list the contents of just ``/data`` after first logging on then your ``/data/te1st`` subdirectory (where ``te1st`` is your username) might not be shown.
However, if you list the contents of ``/data/te1st`` itself or change into that directory
then its contents will appear.  

Later on if you list the contents of ``/data`` again 
you may find that ``/data/te1st`` has disappeared again, as 
it is automatically *unmounted* following a period of inactivity.  

-----

.. _fastdata_dir:

*Fastdata* areas
----------------

**Fastdata** areas are **optimised for large file operations**.  
These areas are `Lustre <https://en.wikipedia.org/wiki/Lustre_(file_system)>`__ filesystems. 

They are are **faster** than :ref:`home_dir`, :ref:`data_dir` and :ref:`shared_dir` when dealing with larger files but 
are **not performant when reading/writing lots of small files** 
(:ref:`scratch_dir` are ideal for reading/writing lots of small temporary files within jobs).
An example of how slow it can be for large numbers of small files is detailed `here <http://www.walkingrandomly.com/?p=6167>`__.

+----------+---------------------+--------+----------------+---------------------+-------------------------+---------------------------+
| System   | Path                | Type   | Quota per user | Filesystem capacity | Shared between systems? | Network bandwith per link |
+==========+=====================+========+================+=====================+=========================+===========================+
| Bessemer | ``/fastdata``       | Lustre | No limits      | 460 TB              | No                      | 25Gb/s Ethernet           |
+----------+---------------------+--------+----------------+---------------------+-------------------------+---------------------------+
| ShARC    | ``/fastdata``       | Lustre | No limits      | 669 TB              | No                      | 100Gb/s (*Omni-Path*)     |
+----------+---------------------+--------+----------------+---------------------+-------------------------+---------------------------+

Snapshotting and mirrored backups
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Snapshotting is not enabled** for fastdata areas and
these areas are **not backed up**.

Managing your files in fastdata areas
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

In order to avoid interference from other users' files 
it is **important** that you store your files in a directory created and named the same as your username 
e.g. if your username is ``te1st`` then you can a fastdata area for yourself using: ::

    mkdir /fastdata/te1st

By default the directory you create will have world-read access.  
If you want to restrict read access to just your account then run ::

    chmod 700 /fastdata/te1st

after creating the directory. 
A more sophisticated sharing scheme would have private and public directories ::

    mkdir /fastdata/te1st
    mkdir /fastdata/te1st/public
    mkdir /fastdata/te1st/private

    chmod 755 /fastdata/te1st
    chmod 755 /fastdata/te1st/public
    chmod 700 /fastdata/te1st/private

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
    **copy important data** from these areas to areas suitable for longer-term storage (:ref:`home_dir`, :ref:`data_dir` (*not* backed up) or :ref:`shared_dir`).

You can use the ``lfs``  command to find out which files in a *fastdata* directory are older than a certain number of days and hence approaching the time of deletion. 
For example, if your username is ``te1st`` then you can find files 50 or more days old using: ::

    lfs find -ctime +50 /fastdata/te1st

File locking
^^^^^^^^^^^^

As of September 2020 POSIX file locking is enabled on all Lustre filesystems. 
Prior to this the lack of file locking support on the University's Lustre filesystems caused problems for certain workflows/applications
(e.g. for programs that create/use SQLite databases).

-----

.. _shared_dir:

*Shared* (project) directories
------------------------------

Each PI at the University is entitled to request a `free 10 TB storage area for sharing data with their group and collaborators <https://sheffield.ac.uk/it-services/research-storage/using-research-storage>`__.
The capacity per area can be extended and additional shared areas can be purchased (both at a cost).

After one of these project storage areas has been requested/purchased it can be accessed in two ways:

* as a Windows-style (SMB) file share on machines other than Bessemer/ShARC using ``\\uosfstore.shef.ac.uk\shared\``;
* as a subdirectory of ``/shared`` on Bessemer/ShARC (you need to **explicitly request HPC access when you order storage from IT Services**).

Snapshotting and mirrored backups
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

+---------------------------+--------------------+----------------------------------------+
| Frequency of snapshotting | Snapshots retained | Backed up onto separate storage system |
+===========================+====================+========================================+
| Every 4 hours             | 10 most recent     | Yes                                    |
+---------------------------+--------------------+----------------------------------------+
| Every night               | Last 7 days        | Yes                                    |
+---------------------------+--------------------+----------------------------------------+

See also: :ref:`recovering_snapshots`.
  
Automounting
^^^^^^^^^^^^

Similar to :ref:`data_dir`, subdirectories beneath ``/shared`` are **mounted on demand** on the HPC systems: 
they may not be visible if you simply list the contents of the ``/shared`` directory but 
will be accessible if you ``cd`` (change directory) into a subdirectory e.g. ``cd /shared/my_group_file_share1``.

Specifics for Bessemer
^^^^^^^^^^^^^^^^^^^^^^

If you need to access a ``/shared`` area on Bessemer please contact `research-it@sheffield.ac.uk <research-it@sheffield.ac.uk>`__ to arrange this.


.. warning::

        * If you access a ``/shared`` directory stored in Sheffield from Bessemer then you may experience slower performance, espeicially for small files.
        * Network traffic between Bessemer and Sheffield Research Filestore is not encrypted when travelling between Sheffield and Leeds over JANET
        * ``/shared`` areas can be created on Bessemer's filestore system if you need faster access from Bessemer

.. _shared_dir_perms:

Permissions behaviour
^^^^^^^^^^^^^^^^^^^^^

You may encounter strange permissions issues when running programs on HPC against the ``/shared`` areas 
e.g. ``chmod +x /shared/mygroup1/myprogram.sh`` fails.
Here we try to explain why.

Behind the scenes, the file server that provides this shared storage manages permissions using 
Windows-style `ACLs <https://en.wikipedia.org/wiki/Access_control_list>`_ 
(which can be set by area owners via the `Research Storage management web interface <https://sheffield.ac.uk/storage/>`__.
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
  a repository is simply moved to ``/shared/someplace`` from elsewhere on Bessemer/ShARC. 
  As a workaround you can tell git to not to track Linux permissions for a single repository using 
  ``git config core.filemode false`` or 
  for all repositories using ``git config --global core.filemode false``.

**Changing how attempts to change permissions are handled**: each ``/shared`` area can be configured so that

#. Attempts to change file/directory mode bits fail (e.g. ``chmod +x /shared/mygroup1/myprogram.sh`` fails) (**default configuration per area**) **or**
#. Attempts to change file/directory mode bits appear to succeed (e.g. ``chmod +x /shared/mygroup1/myprogram.sh`` does not fail but also does not actually change any permissions on the underlying file server) (**alternative configuration per area**)

If you would like to switch to using the second way of handling permissions for a particular ``/shared/`` area then
the Owner of this area should make a request via the Helpdesk.

Further information
^^^^^^^^^^^^^^^^^^^

The documentation for the ``/shared`` storage service includes information on:

* `how access/permissions are managed <https://www.sheffield.ac.uk/it-services/research-storage/access-rights>`__
* `how to create folders with associated permissions <https://www.sheffield.ac.uk/it-services/research-storage/create-folders>`__ 
  within ``/shared`` storage areas

-----

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
The name of this directory is stored in the ``$TMPDIR`` environment variable e.g. on ShARC: ::

    [te1st@sharc-login1 ~]$ qrshx
    [te1st@sharc-node003 ~]$ cd $TMPDIR
    [te1st@sharc-node003 667443.1.all.q]$ pwd
    /scratch/667443.1.all.q

The scheduler will then clean up (delete) ``$TMPDIR`` at the end of your job, 
ensuring that the space can be used by other users.

.. warning::

   If using ``qrsh`` on ShARC to start an interactive job then 
   the ``TMPDIR`` environment variable will unfortunately be undefined
   so you will need to manually create a directory under ``/scratch`` (named using your username)
   and this will not be cleaned up when the job ends.

Anything under the ``/scratch`` may be deleted periodically when the worker-node is idle. 
``/scratch`` is **not backed up**.  There are no quotas for ``/scratch`` storage.

``/scratch`` uses the ext4 filesystem.

-----

.. _community_dir:

*Community* areas for software
------------------------------

Most data that researchers want to share with their collaborators at the University should reside in :ref:`shared_dir`.
However, as mentioned in :ref:`shared_dir_perms`, these areas may not be ideal for storing executable software/scripts
due to the way permissions are handled beneath ``/shared``.

Also, users may want to install software on the clusters that they want to be accessible by all cluster users.

To address these two needs users can request the creation of a new directory beneath of the three directories listed below
and if their request is granted they will be given write access to this area:

+----------+--------------------------+------+-----------------------------+-------------------------------------+-----------------------------------------+
| System   | Path                     | Type | Software install guidelines | Public index of areas               | Notes                                   |
+==========+==========================+======+=============================+=====================================+=========================================+
| Bessemer | ``/usr/local/community`` | NFS  |                             |                                     |                                         |
+----------+--------------------------+------+-----------------------------+-------------------------------------+-----------------------------------------+
| ShARC    | ``/usr/local/community`` | NFS  | :ref:`sharc-community`      | :ref:`sharc-software-install-guide` | Also available at ``/usr/local/extras`` |
+----------+--------------------------+------+-----------------------------+-------------------------------------+-----------------------------------------+

Note that:

* Software installation should follow our installation guidelines where provided.
* Software installations must be maintained by a responsible owner.
* Software which is not actively maintained may be removed.

-----

.. _quota_check:

How to check your quota usage
-----------------------------

To find out your storage quota usage for your :ref:`home directory <home_dir>` and :ref:`data directory <data_dir>` (if not on Bessemer) 
you can use the ``quota`` command:

.. code-block:: console

       [foo11b@sharc-node004 binary]$ quota

       Size  Used Avail Use%  Mounted on
       10G    10G    0G 100%  /home/foo11b
       100G     0  100G   0%  /data/foo11b

In the above, you can see that the quota was set to 10 gigabytes and all of this is in use which is likely to cause jobs to fail.

To determine usage in a particular :ref:`shared_dir` you can use the ``df`` command like so: 

.. code-block:: console

    [foo11b@sharc-node004 binary]$ df -h /shared/myproject1
    Filesystem                        Size  Used Avail Use% Mounted on
    172.X.X.X:/myproject1/myproject1   10T  9.1T  985G  91% /shared/myproject1

To assess what is using up your quota within a given directory, you can make use of the 
:ref:`ncdu module on ShARC <ncdu_sharc>` or the 
:ref:`ncdu module on Bessemer <ncdu_bessemer>`. The **ncdu** utility will give you an 
interactive display of what files/folders are taking up storage in a given directory tree.

-----

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
  or remove them from Bessemer/ShARC completely.

-----

.. _recovering_snapshots:

Recovering files from snapshots
-------------------------------

:ref:`home_dir`, :ref:`data_dir` and :ref:`shared_dir` are regularly :term:`snapshotted <Snapshotted storage>`.
See above for details of the snapshot schedules per area.
A subset of snapshots can be accessed by HPC users from the HPC systems themselves
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
