# LT 1/09/2022

## Cost estimation run on ~5% of partitions
## annotate gnomADv4 with frequencies and variant QC standard fields
## export the rows as a hail.Table on a cpg bucket

import hail as hl
from cpg_utils.hail_batch import output_path

# gnomAD v4 exomes VDS
GNOMAD_V4_0_RAW_EXOMES_VDS_URL = 'gs://gnomad/v4.0/raw/exomes/gnomad_v4.0.vds'

# connect to hail, requester_pays parameters are setup in the dataproc parameters
hl.init(default_reference='GRCh38')

# read vds
big_vds = hl.vds.read_vds(GNOMAD_V4_0_RAW_EXOMES_VDS_URL)

# subset to ~5% of partitions (big_vds has 47960 partitions)
small_vds = hl.vds.VariantDataset(
    big_vds.reference_data._filter_partitions(range(2500)),
    big_vds.variant_data._filter_partitions(range(2500)),
)

# densify as a MatrixTable
mt = hl.vds.to_dense_mt(small_vds)

# remove REF/REF rows left after densification
mt = mt.filter_rows(hl.len(mt.alleles) < 2, keep=False)

# get alleles frequencies
mt = mt.transmute_entries(GT=hl.experimental.lgt_to_gt(mt.LGT, mt.LA))
mt = hl.variant_qc(mt)

# keep only the rows and repartition
# given the drastic size reduction by summarizing the columns, reduce the number of partitions ~10 fold
ht = mt.rows().repartition(250, shuffle=False)

# write table to cpg bucket
ht.write(output_path('gnomad_v4.0_raw_exomes.costrun.ht'), overwrite=True)
