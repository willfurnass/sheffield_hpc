.. _troubleshooting:

Troubleshooting
===============
In this section, we'll discuss some tips for solving problems with iceberg. It is suggested that you work through some of the ideas here before contacting the service desk for assistance.

Exceeding your disk space quota
-------------------------------
Each user of the system has a fixed amount of disk space available in their home directory. If you exceed this quota, various problems can emerge such as an inability to launch applications and run jobs.
To see if you have exceeded your disk space quota, run the quota command:

.. code-block:: none

       quota

       Size  Used Avail Use% Mounted on
       5.1G  5.1G     0 100% /home/foo11b
        50G     0   50G   0% /data/foo11b

In the above, you can see that the quota was set to 5.1 gigabytes and all of this is in use. Any jobs submitted by this user will likely result in an Eqw status. The recommended action is for the user to delete enough files to allow normal work to continue.
