# make_transcript_coord_file.R
# M. Silk


# Given an ENST transcript, generate
# GRCh37/hg19 coordinates for use in
# the structural constraint prototyping


library(data.table)
library(magrittr)
library(purrr)
library(rtracklayer)
library(seqinr)


# Read in GTF file, filter to
read_gtf <- function(fname) {
  rtracklayer::import(fname) %>%
    .[mcols(.)$type %in% "CDS"]
}


# Make coordinates sequence table
make_sequence_table <- function(my_enst) {

  get_cds_ranges <- function(my_enst) {
    gtf[!is.na(mcols(gtf)$transcript_id)] %>%
      .[startsWith(mcols(.)$transcript_id, my_enst)]
  }

  get_coding_sequence <- function(my_enst) {
    enst_coding[enst == my_enst]$cds
  }

  remove_stop_codon <- function(s) {
    substr(s, 1, nchar(s) - 3)
  }

  expand_range <- function(start, end) {
    start:end
  }

  get_chromosome <- function(my_gtf) {
    chrom(my_gtf) %>%
      as.character() %>%
      unique()
  }

  get_strand <- function(my_gtf) {
    strand(my_gtf) %>%
      as.character() %>%
      unique()
  }

  orient_coding_sequence <- function(my_gtf, my_coding) {
    if (get_strand(my_gtf) == "+") {
      return(my_coding %>% strsplit("") %>% .[[1]])
    } else {
      return(my_coding %>% strsplit("") %>% .[[1]] %>% rev() %>% map(invert_base) %>% unlist())
    }
  }

  orient_protein_pos <- function(my_gtf, my_coding) {
    if (get_strand(my_gtf) == "+") {
      return(seq_len(nchar(my_coding) / 3) %>% rep(each = 3))
    } else {
      return(seq_len(nchar(my_coding) / 3) %>% rep(each = 3) %>% rev())
    }
  }

  cds_ranges <- get_cds_ranges(my_enst)
  cds <- get_coding_sequence(my_enst) %>%
    remove_stop_codon()

  stopifnot(identical(
    width(cds_ranges) %>% sum(),
    nchar(cds)
  ))

  data.table(
    chrom = get_chromosome(cds_ranges),
    pos = get_positions(cds_ranges),
    ref = orient_coding_sequence(cds_ranges, cds),
    aapos = orient_protein_pos(cds_ranges, cds)
  )
}


run_transcript <- function(my_enst) {

  x = make_sequence_table(my_enst)

  # Add additional columns
  extra_cols <- idmapping[enst == my_enst, .(uniprot_id = uniprot,
                                             hgnc = gene,
                                             ensg,
                                             transcript_id = enst,
                                             assembly = "GRCh37")]

  x = cbind(x, extra_cols)

  # Create reference position column
  x[, genomic_position := paste0(chrom %>% gsub("chr", "", .),
                                 ":",
                                 pos)]

  # Reformat
  x <- x[, .(uniprot_id,
             hgnc,
             ensg,
             transcript_id,
             protein_position = aapos,
             assembly,
             genomic_position,
             ref)]

  # Write to file
  fwrite(x, paste0("alignment/positions/positions_", x$hgnc[1], "_", my_enst, ".txt.gz"))
}


# Read in data
idmapping <- fread("alignment/gencode_translations_matched.txt.gz")
gtf <- read_gtf("alignment/gencode.v19.annotation.gtf.gz")
enst_coding <- fread("alignment/gencode_transcripts.txt.gz")


# Run for all transcripts with a matching ENST to AF2 structure
valid_transcripts <- idmapping[match == TRUE]$enst %>% unique
map(valid_transcripts[1:10], run_transcript)
