# LT 31/08/2022

## annotate gnomADv4 with frequencies and variant QC standard fields
## export the rows as a hail.Table on a cpg bucket

import hail as hl
from cpg_utils.hail_batch import output_path

# connect to hail, requester_pays parameters are setup in the dataproc parameters
hl.init(default_reference='GRCh38')

# gnomAD v4 exomes VDS
vds = hl.vds.read_vds('gs://gnomad/v4.0/raw/exomes/gnomad_v4.0.vds')

# densify as a MatrixTable
mt = hl.vds.to_dense_mt(vds)

# remove REF/REF rows left after densification
mt = mt.filter_rows(hl.len(mt.alleles) < 2, keep=False)

# get alleles frequencies
mt = mt.transmute_entries(GT=hl.experimental.lgt_to_gt(mt.LGT, mt.LA))
mt = hl.variant_qc(mt)

# keep only the rows and repartition
# given the drastic size reduction by summarizing the columns, reduce the number of partitions ~50 fold
ht = mt.rows().naive_coalesce(1000)

# write table to cpg bucket
ht.write(output_path('gnomad_v4.0_raw_exomes.ht'), overwrite=True)
