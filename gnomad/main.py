#!/usr/bin/env python3

"""First run of a dataproc in the USA for gnomADv4 analysis"""

import hailtop.batch as hb
from cpg_utils.hail_batch import get_config, remote_tmpdir, output_path
from analysis_runner import dataproc
import subprocess

config = get_config()

service_backend = hb.ServiceBackend(
    billing_project=config['hail']['billing_project'], remote_tmpdir=remote_tmpdir()
)

batch = hb.Batch(name='dataproc gnomADv4', backend=service_backend)


cluster = dataproc.setup_dataproc(
    batch,
    max_age='100h',
    #packages=['click', 'selenium'], # is this needed ? NO
    #packages=['çpg-utils']
    init=['gs://cpg-reference/hail_dataproc/install_common.sh'], # what does this do ? plotting, analysis-runner
    cluster_name='gnomADv4',
    worker_machine_type = 'n1-standard-8',
    num_workers = 20,
    worker_boot_disk_size = 50,
    num_secondary_workers = 0,
    #secondary_worker_boot_disk_size = 20, # is this needed ?
    region='us-central1',
    requester_pays_allow_all=True
)

cluster.add_job('gnomadv4_extract_test.py', job_name='first_v4_extract_test')

# add the needed files
subprocess.run(['gsutil', 'cp', 'transcripts.json', output_path('transcripts.json')], check=True)


# Don't wait, which avoids resubmissions if this job gets preempted.
batch.run(wait=False)