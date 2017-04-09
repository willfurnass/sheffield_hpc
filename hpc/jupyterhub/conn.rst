.. _jh_conn: 

JupyterHub on ShARC: Connecting
===============================

#. From a web browser navigate to the following URL 
   from a machine that is on-campus or
   has a `VPN connection to the University <https://www.sheffield.ac.uk/cics/vpn>`_::

        https://jupyter-sharc.shef.ac.uk

#. When prompted, log in using your ShARC/Iceberg username 
   (e.g. ``te1st``, *not* the first part of your email address)
   and password.

    .. image:: /images/jupyterhub/jh-sharc-login.png

#. Click **Start server** when prompted:

    .. image:: /images/jupyterhub/jh-sharc-start-server.png

#. Wait for your Jupyter session to start on one of the cluster's worker nodes:

    .. image:: /images/jupyterhub/jh-sharc-server-starting.png

.. warning::
    If the cluster is busy **or** 
    you have requested an incompatible or otherwise unsatisfiable set of resources 
    from the job scheduler
    then this attempt to start a session will time out
    and you will return to... **CHECK EXACTLY WHAT HAPPENS**


**TODO: FINISH - CAN ONLY START ONE SESSION - INC REF TO NEXT SECTION**

MAIN INTERFACE: 

.. image:: /images/jupyterhub/sharc-jh-main-nb-svr-interface.png
