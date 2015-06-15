Modules
-------

In general the software available on iceberg is loaded and unloaded via the use
of the modules system [#env-modules]_.


Available modules can be listed using the::

    module avail

command, modules can be loaded by running::

    module load apps/python/2.7

and unloaded using::

    module unload apps/python/2.7

It is possible to load multiple modules at once, to create your own environment
with just the software you need. For more details on how modules work run the 
`man module` command at an iceberg command prompt.


.. [#env-modules] http://modules.sourceforge.net/
