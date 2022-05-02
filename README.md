# structural-constraint
Calculating missense constraint within protein tertiary structures. The following provides information on how to run each script.

## gencode_translations_matched.txt.gz - Additional files

You will also need to download the following, and set the paths within the R script accordingly:

### LiftOver executable
You will need to make a free educational account and sign in to download the executable.
Add to cart the LiftOver program https://genome-store.ucsc.edu, click Login/Register in the top right, follow the instructions then check out.

### LiftOver chain file
`wget -c http://hgdownload.cse.ucsc.edu/goldenPath/hg38/liftOver/hg38ToHg19.over.chain.gz`

### GENCODE GTF file
`wget -c http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_39/gencode.v39.annotation.gtf.gz`
