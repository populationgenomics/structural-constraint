#!/usr/bin/env python3

"""COST RUN: Extract rows fron gnomAD v4"""

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
    max_age='30h', # increased for safety given decrease of #workers
    packages=['cpg-utils'],
    init=[],
    cluster_name='gnomADv4',
    worker_machine_type = 'n1-standard-8',
    num_workers = 2,
    # note: the image requires at least 30 GB on boot disk
    worker_boot_disk_size = 50, # total HDFS space 100 GB, HDFS peak-use during the cost test run was 21 GiB
    num_secondary_workers = 10, # total vCPUs ~ 80 + 16, decreased compared to cost test run to reduce inactive secondary workers during table write
    region='us-central1',
    requester_pays_allow_all=True
)

cluster.add_job('gnomadv4_cost_run_extract_rows.py', job_name='gnomadv4_cost_run_extract_rows')

# Don't wait, which avoids resubmissions if this job gets preempted.
batch.run(wait=False)
