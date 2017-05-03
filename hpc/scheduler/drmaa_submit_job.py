#!/usr/bin/env python

from __future__ import print_function
import drmaa
import os


def main():
    """
    Create a DRMAA session then submit a job.
    Note, need file called myjob.sh in the current directory.
    """
    with drmaa.Session() as s:
        print('Creating a job template')
        jt = s.createJobTemplate()
        # The job is to run an executable in the current working directory
        jt.remoteCommand = os.path.join(os.getcwd(), 'myjob.sh')
        # Arguments to the remote command
        jt.args = ['alpha', 'beta']
        # Join the standard output and error logs
        jt.joinFiles = True

        job_id = s.runJob(jt)
        print('Job {} submitted'.format(job_id))

        print('Cleaning up')
        s.deleteJobTemplate(jt)


if __name__ == '__main__':
    main()
