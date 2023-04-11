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
    .[mcols(.)$type %in% c("CDS")]
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

  get_positions <- function(my_gtf) {
    map2(start(my_gtf), end(my_gtf), expand_range) %>%
      unlist() %>%
      sort()
  }

  invert_base <- function(base) {
    if (base == "A") {
      return("T")
    } else if (base == "C") {
      return("G")
    } else if (base == "G") {
      return("C")
    } else {
      return("A")
    }
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
  cds <- get_coding_sequence(my_enst)
  cds_ranges_width <- width(cds_ranges) %>% sum()
  cds_width <- nchar(cds)

  # Remove stop codons where length of the coding sequence
  # is 3 less than the ranges (almost all except oddballs)
  if (cds_ranges_width == cds_width - 3) {
    cds <- cds %>% remove_stop_codon
  }

  # Stop if coding sequence not divisible by 3 (oddballs)
  stopifnot((width(cds_ranges) %>% sum) %% 3 == 0)

  # Stop if lengths don't match
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
  fwrite(x, paste0("alignment/positions_canonical/positions_", x$uniprot_id[1], "_", x$hgnc[1], "_", my_enst, ".tsv.gz"), sep = "\t")
}


# Read in data
idmapping <- fread("alignment/gencode_translations_matched.txt.gz")
gtf <- read_gtf("alignment/gencode.v19.annotation.gtf.gz")
enst_coding <- fread("alignment/gencode_transcripts.txt.gz")


# Run for all transcripts with a matching ENST to AF2 structure
valid_transcripts <- idmapping[match == TRUE]$enst %>% unique


# Canonical transcripts: any with "basic" as the designation
# in the GENCODE GTF files. Strip the transcript version at the end
# and filter valid_transcripts by this set
gtf_basic_transcripts = gtf[which(gtf$tag == "basic")]$transcript_id %>%
  unique %>%
  strsplit(".", fixed = TRUE) %>%
  map(1) %>%
  unlist
canonical_transcripts = valid_transcripts[valid_transcripts %in% gtf_basic_transcripts]


# Wrapper function to allow error messages saved as a side-effect,
# similar to Haskell's Maybe. Returns the results of run_transcript
# in the first item of each list element, error messages in the second.
safely_run_transcript <- safely(.f = run_transcript)
map(valid_transcripts, safely_run_transcript)
a = map(canonical_transcripts, safely_run_transcript)

# there are about 150 transcripts that don't work which we will leave for now
# less issues with the newest set though





# more testing
gtf_transcripts = gtf$transcript_id %>% unique %>% strsplit(".", fixed = TRUE) %>% map(1) %>% unlist

# how many in idmapping
gtf_transcripts[gtf_transcripts %in% idmapping$enst] %>% length

# almost all of the GTF basic transcripts are in idmapping (15,080)
# but then none are in valid transcripts?
idmapping[match == TRUE & enst %in% canonical_transcripts]
idmapping[enst %in% canonical_transcripts]
idmapping[enst %in% gtf_basic_transcripts]


# appris
# note: no version numbers to deal with
appris = fread("~/data/appris/appris_data.principal_gencode19ensembl74.txt", header = FALSE)
idmapping[match == TRUE & enst %in% appris[V5 == "PRINCIPAL:1"]$V3]
# and that completely restores the set of transcripts
idmapping[match == TRUE & enst %in% appris[V5 == "PRINCIPAL:1"]$V3]$gene %>% unique %>% length
# 13,943 matching genes

# compare with mane
mane = fread("~/Desktop/MANE.GRCh38.v0.95.summary.txt")
mane_transcripts = mane$Ensembl_nuc %>% strsplit(".", fixed = TRUE) %>% map(1) %>% unlist %>% unique

gtf_basic_transcripts[gtf_basic_transcripts %in% mane_transcripts]




# delete later
x = rtracklayer::import("~/data/gencode/gencode.v39.annotation.gtf")
