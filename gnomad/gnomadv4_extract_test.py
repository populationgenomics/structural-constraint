# LT 23/08/2022

## Extract a small subset of 10 transcripts from gnomADv4
## annotate with frequencies
## export the rows as a hail.Table


# target transcripts for panel of initial 10 genes selected by the structural constraint project
# read from file 'transcripts.json'

import json
import hail as hl
from cpg_utils.hail_batch import output_path

transcripts_file = open("transcripts.json", "r")
TRANSCRIPTS = json.load(transcripts_file)

# connect to hail, using appropriate requester_pays setup

hl.init(default_reference='GRCh38',
        spark_conf = {
            'spark.hadoop.fs.gs.requester.pays.mode': 'CUSTOM',
            'spark.hadoop.fs.gs.requester.pays.buckets': 'gnomad',
            'spark.hadoop.fs.gs.requester.pays.project.id': 'constraint-232598',
        })

# gnomAD v4 exomes VDS
big_vds = hl.vds.read_vds('gs://gnomad/v4.0/raw/exomes/gnomad_v4.0.vds')

# find genomics coordinates of target transcripts and filter vds
# use GENCODE v39 to match Alphafold 2


# to generate a blocked zipped file I had to
#
# curl https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_39/gencode.v39.basic.annotation.gtf.gz --output ./gencode.v39.basic.annotation.gtf.gz
# sudo apt update
# sudo apt install tabix
# zcat gencode.v39.basic.annotation.gtf.gz | bgzip -c > gencode.v39.basic.annotation.gtf.bgz
GENCODE = 'gencode.v39.basic.annotation.gtf.bgz'
target_intervals = hl.experimental.get_gene_intervals(transcript_ids=TRANSCRIPTS, reference_genome='GRCh38', gtf_file=GENCODE)
small_vds = hl.vds.filter_intervals(big_vds, target_intervals, split_reference_blocks=False, keep=True)

# densify as a MatrixTable
mt = hl.vds.to_dense_mt(small_vds)

# get alleles frequencies
mt = mt.transmute_entries(GT=hl.experimental.lgt_to_gt(mt.LGT, mt.LA))
mt = hl.variant_qc(mt)

# keep only the rows
ht = mt.rows()

# safe guard, before writting to an Australian bucket. Limit to 1 000 000 variants
# note that there is no need to write to Australia, since this is an intermediary result
# it would be preferable to have a bucket hosted in the USA
if (ht.count() < 1e6) :
    # use cpg utils to find proper path (bucket will change with access level)
    ht.write(output_path("small_ht"), overwrite=True)


