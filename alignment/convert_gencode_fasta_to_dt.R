
# convert_gencode_fasta_to_dt.R
# M. Silk

# Convert GENCODE FASTA to data.table
# FASTA fields in the header are variable in
# what's included (in content, order and length)


library(data.table)
library(magrittr)
library(purrr)
library(seqinr)


gencode_translations <- read.fasta("alignment/gencode.v19.pc_translations.fa.gz")


# If running on modern GENCODE pc_translations files,
# colnames must also include "ensp" as the first field
# This is NOT included in GENCODE v19
colnames_translations <- c(
  "enst", "ensg",
  "havana_g", "havana_t", "transcript",
  "gene", "aalength", "sequence"
)


# Convert FASTA entry into header columns attached to
# sequence as a single string
convert_entry <- function(x) {
  header <- seqinr::getName(x) %>%
    strsplit("|", fixed = TRUE) %>%
    .[[1]]
  s <- seqinr::getSequence(x) %>%
    unlist() %>%
    as.character() %>%
    toupper() %>%
    paste(collapse = "")
  vals <- c(header, s) %>%
    as.list() %>%
    as.data.table()
  names(vals) <- colnames_translations
  return(vals)
}


gencode_translations_dt <- map(gencode_translations, convert_entry) %>%
  rbindlist(use.names = TRUE, fill = TRUE)

fwrite(gencode_translations_dt, "alignment/gencode_translations.txt.gz", sep = "\t")
