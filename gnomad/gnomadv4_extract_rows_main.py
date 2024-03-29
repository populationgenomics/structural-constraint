#!/usr/bin/env python3

"""Extract rows from gnomAD v4"""

import hailtop.batch as hb
from cpg_utils.hail_batch import get_config, remote_tmpdir
from analysis_runner import dataproc
import subprocess

config = get_config()

service_backend = hb.ServiceBackend(
    billing_project=config['hail']['billing_project'], remote_tmpdir=remote_tmpdir()
)

batch = hb.Batch(name='dataproc gnomADv4', backend=service_backend)


cluster = dataproc.setup_dataproc(
    batch,
    max_age='330h', # increased for safety
    packages=['cpg-utils'],
    init=[],
    cluster_name='gnomADv4',
    worker_machine_type = 'n1-highmem-8',
    num_workers = 2,
    master_boot_disk_size = 1000, # increase disk space even more, as a potential fix for dataproc crash
    worker_boot_disk_size = 1000, # total HDFS space 2TB, HDFS peak use during test run was 450 GiB
    secondary_worker_boot_disk_size = 500, # increase disk space on secondary workers
    num_secondary_workers = 100, # total vCPUs ~ 800, for a result Table with 1000 partitions
    autoscaling_policy = 'max100', # allows decommission of useless workers
    region='us-central1',
    requester_pays_allow_all=True,
    stop_cluster=False # for debugging cluster crash
)

cluster.add_job('gnomadv4_extract_rows.py', job_name='gnomadv4_extract_rows')

# Don't wait, which avoids resubmissions if this job gets preempted.
batch.run(wait=False)