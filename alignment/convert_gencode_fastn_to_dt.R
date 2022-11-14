
# convert_gencode_fasta_to_dt.R
# M. Silk

# Convert GENCODE FASTA to data.table
# FASTA fields in the header are variable in
# what's included (in content, order and length)


library(data.table)
library(magrittr)
library(purrr)
library(seqinr)


gencode_transcripts <- read.fasta("alignment/gencode.v19.pc_transcripts.fa.gz")

colnames_transcripts = c("enst", "ensg",
                         "havana_g", "havana_t", "transcript",
                         "gene", "ntlength",
                         "index1", "index2", "index3",
                         "sequence")


# Column numbers are variable here. Some have
# 5' and 3' UTR regions to the transcripts, frustratingly
convert_entry <- function(a) {
  header = attr(a, "name") %>% strsplit("|", fixed = T) %>% .[[1]]
  s = getSequence(a) %>% unlist %>% toupper %>% paste(collapse = "")
  vals = c(s, header)
  vals %<>% as.list %>% as.data.table
  return(vals)
}


gencode_transcripts_dt = lapply(gencode_transcripts, convert_entry) %>%
  rbindlist(use.names = T, fill = T)


# Reorder columns
gencode_transcripts_dt %<>% .[, c(2:11, 1)]


# Add names
names(gencode_transcripts_dt) = colnames_transcripts


# Trim transcript versions from ID's
gencode_transcripts_dt[, enst := strsplit(enst, ".", fixed = TRUE)[[1]][1],
                       seq_len(nrow(gencode_transcripts_dt))]


# Trim sequences to just CDS range
# Target index may be in the 1st or 2nd column,
# select whichever starts with "CDS"
trim_seqs_to_cds <- function(my_seqs) {

  # Edit index string formatting from seqs to
  # get just CDS ranges as vector,
  # eg. from "CDS:104-3678" --> 104 3678
  edit_index_string <- function(s) {
    substr(s, 5, nchar(s)) %>%
      strsplit("-") %>%
      .[[1]] %>%
      as.numeric()
  }

  # Trim sequences to only the CDS range
  # Given a sequence s, trim it using
  # index string of ranges
  trim_sequence <- function(s, i) {
    range <- edit_index_string(i)
    return(substr(s, range[1], range[2]))
  }

  my_seqs[startsWith(index1, "CDS"), cds := map2(sequence, index1, trim_sequence) %>% unlist()]
  my_seqs[startsWith(index2, "CDS"), cds := map2(sequence, index2, trim_sequence) %>% unlist()]
}


gencode_transcripts_dt <- trim_seqs_to_cds(gencode_transcripts_dt)


# Remove index columns and original sequence
gencode_transcripts_dt[, sequence := NULL]
gencode_transcripts_dt[, index1 := NULL]
gencode_transcripts_dt[, index2 := NULL]
gencode_transcripts_dt[, index3 := NULL]


# Write to file
fwrite(gencode_transcripts_dt, "alignment/gencode_transcripts.txt.gz")
