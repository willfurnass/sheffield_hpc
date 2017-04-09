.. _jh_issues: 

JupyterHub: errors and troubleshooting
======================================









NOTES BELOW
-----------

Less common issues:

#. If you modify the ``PYTHON_PATH`` variable in your ``.bashrc`` file your
   Jupyter server may not start correctly, and you may encounter a 503 error
   after logging into the hub. The solution to this is to remove these lines
   from your ``.bashrc`` file.
#. If you have previously tried installing and running Jupyter yourself (i.e.
   not using this JupyterHub interface) then you may get 503 errors when
   connecting to JupyterHub due to `the old .jupyter profile in your home
   directory <https://github.com/jupyter/jupyterhub/issues/294>`_;  if you then
   find a Jupyter log file in your home directory containing ``SSL
   WRONG_VERSION_NUMBER`` then try deleting (or renaming) the ``.jupyter``
   directory in your home directory.

#. CAN'T DELETE ENVS?
#. CAN'T BROWSE OUT OF HOME DIR?
