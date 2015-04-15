iceberg Documentation
=====================

This repository builds the iceberg documentation using sphinx.

Currently, it harvests the software documentation from
`/usr/local/packages6` and `/usr/local/extras` so needs to be run on
a worker node on iceberg. 

You can activate the python environment to use sphinx by running:

```
source activate_env.sh
```

once this has been activated you can run 

```
make html
```

to build the documentation into a website, you can then view the
documentation by running:

```
firefox _build/html/index.html
```

because sphinx puts the output files into the `_build/html`
directory.


Deploying
---------

The files in `_build/html` can be copied to Cpanel to deploy the
site.
