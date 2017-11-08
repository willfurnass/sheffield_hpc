.. _jh_sess_monit: 

Monitoring and controlling your Jupyter session
===============================================

If you want to see what Jupyter Notebooks and Terminals your Jupyter server is running then

#. Click the *Home* browser tab;
#. Click the *Running* browser tab.

Here you can **shut down (close) individual Notebooks**, which may be useful to free up 
the memory and CPU cores that the cluster's job scheduler has alllocated for your JupyterHub session.

.. image:: /images/jupyterhub/sharc-jh-show-running.png

----

If you want to **stop your entire Jupyter server** (i.e. end your JupyterHub session) you can 

#. Click *Control Panel* on the *Home* browser tab;
#. Click *Stop My Server*

----

Your JupyterHub session **may terminate for other reasons**:

* The cluster's job scheduler may stop your JupyterHub job if your Jupyter server exceeds the amount of RAM you requested via the :ref:`Spawner Options <jh_spawner_opts>` page.
* Your JupyterHub job has been running for longer than the (fixed) duration specified on the :ref:`Spawner Options <jh_spawner_opts>` page.

