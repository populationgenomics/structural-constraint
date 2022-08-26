# LT 23/08/2022

## Extract a small subset of 10 transcripts from gnomADv4
## annotate with frequencies
## export the rows as a hail.Table


# target transcripts for panel of initial 10 genes selected by the structural constraint project
# read from file 'transcripts.json'

import json

transcripts_file = open("transcripts.json", "r")
TRANSCRIPTS = json.load(transcripts_file)

# connect to hail, using appropriate requester_pays setup
import hail as hl
hl.init(default_reference='GRCh38',
        spark_conf = {
            'spark.hadoop.fs.gs.requester.pays.mode': 'CUSTOM',
            'spark.hadoop.fs.gs.requester.pays.buckets': 'gnomad',
            'spark.hadoop.fs.gs.requester.pays.project.id': 'constraint-232598',
        })

# gnomAD v4 exomes VDS
big_vds = hl.vds.read_vds('gs://gnomad/v4.0/raw/exomes/gnomad_v4.0.vds')

# find genomics coordinates of target transcripts and filter vds
target_intervals = hl.experimental.get_gene_intervals(transcript_ids=TRANSCRIPTS, reference_genome='GRCh38')
small_vds = hl.vds.filter_intervals(big_vds, target_intervals, split_reference_blocks=False, keep=True)

# densify as a MatrixTable
mt = hl.vds.to_dense_mt(small_vds)

# get alleles frequencies
mt = mt.transmute_entries(GT=hl.experimental.lgt_to_gt(mt.LGT, mt.LA))
mt = hl.variant_qc(mt)

# keep only the rows and write down
ht = mt.rows()
ht.write("small_ht", overwrite=True)

