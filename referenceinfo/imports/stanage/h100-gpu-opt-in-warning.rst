.. warning::

    During the 2 week introduction phase of the H100 GPUs to the Stanage cluster, usage of the H100 GPUs requires the ``--partition=gpu-h100`` and ``--gres=gpu:1`` arguments to be set in your submission scripts.
    This is to ensure usage is "opt in" by users as the slightly different architecture of these GPUs to the existing A100 GPUs may necessitate changes to batch submission scripts and selected software versions.

    Eventually the H100 GPUs will be brought into the general GPU partition, at which point the ``--partition=gpu`` will be required to access H100s (or any GPUs). At that stage any submissions
    using the general ``--gres=gpu:1`` will be scheduled with the first available GPU of any type. Requesting a specific type of GPU will then require selection via the ``--gres=gpu:h100:1`` or
    ``--gres=gpu:a100:1`` arguments.

    When these latter changes are made, we will give advanced notice via email and by amendments made within this documentation.
