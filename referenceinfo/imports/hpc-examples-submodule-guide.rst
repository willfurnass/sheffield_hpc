How to Make Changes to ``hpc-examples`` and Reflect Them in ``sheffield_hpc``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

If you make changes inside the ``hpc-examples`` submodule (e.g. editing or adding a script), follow these steps to ensure both repositories stay in sync:

First, commit and push the changes in ``hpc-examples``:

.. code-block:: bash

   cd hpc-examples
   git add <file>
   git commit -m "Add/update example script"
   git push

Wait for the ``hpc-examples`` pull request to be reviewed, approved, and merged into ``main`` before updating the submodule pointer in ``sheffield_hpc``:

.. code-block:: bash

   cd ..
   git add hpc-examples
   git commit -m "Update hpc-examples submodule pointer"
   git push

This ensures ``sheffield_hpc`` tracks a public, stable commit of ``hpc-examples``.

.. important::

   When updating the submodule pointer to ``hpc-examples``:

   - Only reference commits that have already been merged into the ``main`` branch.
   - Never commit a pointer to a branch-only or unpublished commit — this will cause submodule fetch failures in GitHub Actions workflows and for collaborators cloning the repo.
   - Local changes to ``hpc-examples`` will build and render correctly within your local clone of ``sheffield_hpc``, even before they are pushed or merged.
   - You may temporarily point to a local commit of ``hpc-examples`` for testing, **but do not push that pointer** to ``sheffield_hpc`` until the referenced commit is available in the upstream ``hpc-examples/main`` branch.
