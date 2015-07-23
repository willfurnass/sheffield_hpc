Sun Grid Engine
===============


Frequently Asked SGE questions
------------------------------
**How do you ensure that a job starts after a specified time?**

Add the following line in your submission script ::

    #$ -a time

but replace ``time`` with a time in the format MMDDhhmm

For example, for 22nd July at 14:10, you’d do ::

    #$ -a 07221410

This won’t guarantee that it will run precisely at this time since that depends on available resources. It will, however, ensure that the job runs *after* this time. If your resource requirements aren’t too heavy, it will be pretty soon after. When I tried it, it started about 10 seconds afterwards but this will vary.
