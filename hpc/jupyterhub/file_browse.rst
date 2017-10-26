.. _jh_file_browse: 

Jupyter's file browser
======================

After starting a Jupyter session you are presented with Jupyter's file browser view,
which is the default tab in Jupyter's user interface.  
This view shows you **files on the machine running your Jupyter session** (here, the cluster), 
*not* your local machine.  This behaves much like a desktop file browser application:

    .. image:: /images/jupyterhub/sharc-jh-main-nb-svr-interface.png

* Click on a directory to browse into it.
* Click on a text file (e.g. a R or Python script) to view/edit it in your browser.  
  You'll notice that Jupyter's text editor understands different programming languages and
  will use colour to make key syntactic elements easier to identify.
* Click on a Notebook (``.ipynb`` filename suffix) to open it in a new browser tab.
* Rename, delete or move files/folders by checking the relevant tick-boxes to the left of each item of interest then
  click the appropriate option above.
* The drop-down menu at the top of the column of check-boxes allows you to filter to view 
  just Notebooks, 
  all Notebooks that are currently running, 
  all files or 
  all folders.

.. _jh_automount_issue:

.. warning:: 

   Certain directories are may not be accessible via this interface:

   * ``/home/username``
   * ``/data/username``
   * ``/shared/volname``

   This set of directories are :ref:`automounted <filestore>` 
   i.e. made available to the user on demand
   but you cannot express that demand via this interface.
   If you browse into ``/data`` and it is empty or does not contain your personal subdirectory 
   then you need to briefly open a :ref:`Jupyter terminal <jh_terminal>` and 
   run: ::

      ls /data/username

   then that directory should subsequently be visible/accessible in this file browser.
