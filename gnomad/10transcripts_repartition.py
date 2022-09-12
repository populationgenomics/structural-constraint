# LT 12/09/2022

## repartition 10 transcripts test dataset from gnomADv4

import hail as hl
from cpg_utils.hail_batch import output_path

hl.init(default_reference='GRCh38')

# gnomAD v4 10 transcripts test table
ht = hl.read_table('gs://cpg-constraint-main/gv4_10transcripts_6/small.ht/')

# repartition and save in test
ht.repartition(100, shuffle=True).write('gs://cpg-constraint-test/gv4_10transcripts.ht')
