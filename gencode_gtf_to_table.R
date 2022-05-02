# gencode_gtf_to_table.R
# M. Silk

# Get ranges for a given transcript based on a gene symbol ID
# Target gene examples ACTA1 LRRK2 FBXW7
#
# GTF data contains one row per datatype, described in the "type"
# column. For instance filtering to an ENST will give you a number
# of rows with type exon, cds, start_codon, transcript, etc etc.
# CDS ranges do not include the stop codon, however the coding
# sequences do, so the stop_codon type is also captured.
#
# idmapping describes the transcripts that exact match the AF2 data
# This has been generated separately, and requires the full AF2 datset.
#
# enst_coding contains strings of the sequence for each transcript
# based on the transcripts relevant to the AF2 dataset.
#
# LiftOver required to get hg19 coords for GRCh37 for now.

library(data.table)
library(magrittr)
library(purrr)
library(rtracklayer)
library(seqinr)

# ===============
# CONSTANTS
GTF_TYPE_FILTER = c("CDS")
GENOME_BUILD = "GRCh37"

LIFTOVER_EXECUTABLE = "/home/msilk/software/liftOver/liftOver"
LIFTOVER_CHAINFILE = "/home/msilk/software/liftOver/hg38ToHg19.over.chain.gz"
UNIPROT_GENCODE_IDMAPPINGS = "/home/msilk/scripts/alphavis/gencode_translations_matched.txt"  # **
UNIPROT_IDS = "/home/msilk/scripts/alphavis/uniprot_enst_matches.txt"                         # **
GENCODE_GTF = "/home/msilk/data/gencode/gencode.v39.annotation.gtf"
GENCODE_TRANSCRIPTS = "/home/msilk/data/gencode/gencode.v39.pc_transcripts.fa"
ENST_CODING = "/home/msilk/scripts/alphavis/enst_coding.txt"                                  # **

# ** custom files from AlphaFold transcript mapping project

# ===============


# ===============
# READ IN DATA
read_idmapping <- function(fname) {
  fread(fname)
}


read_uniprot_ids <- function(fname) {
  fread(fname)
}


read_gtf <- function(fname) {
  rtracklayer::import(fname) %>%
    .[mcols(.)$type %in% GTF_TYPE_FILTER]
}


read_coding_translations <- function(fname) {
  fread(ENST_CODING)
}


get_enst <- function(my_gene) {
  idmapping[gene == my_gene & matchcheck == "Match" & mane_select == TRUE]$enst
}


get_cds_ranges <- function(my_enst) {
  gtf[!is.na(mcols(gtf)$transcript_id)] %>%
    .[startsWith(mcols(.)$transcript_id, my_enst)]
}


get_coding_sequence <- function(my_enst) {
  enst_coding[ensembl_transcript_id == my_enst]$coding
}


remove_stop_codon <- function(s) {
  substr(s, 1, nchar(s) - 3)
}


expand_range <- function(start, end) {
  start:end
}


get_positions <- function(my_gtf) {
  map2(start(my_gtf), end(my_gtf), expand_range) %>% unlist %>% sort
}


get_chromosome <- function(my_gtf) {
  chrom(my_gtf) %>% as.character %>% unique
}


get_strand <- function(my_gtf) {
  strand(my_gtf) %>% as.character %>% unique
}


invert_base <- function(base) {
  if (base == "A") {
    return("T")
  }
  else if (base == "C") {
    return("G")
  }
  else if (base == "G") {
    return("C")
  }
  else {
    return("A")
  }
}

orient_coding_sequence <- function(my_gtf, my_coding) {
  if (get_strand(my_gtf) == "+") {
    return(my_coding %>% strsplit("") %>% .[[1]])
  }
  else {
    return(my_coding %>% strsplit("") %>% .[[1]] %>% rev %>% map(invert_base) %>% unlist)
  }
}


orient_protein_pos <- function(my_gtf, my_coding) {
  if (get_strand(my_gtf) == "+") {
    return(seq_len(nchar(my_coding) / 3) %>% rep(each = 3))
  }
  else {
    return(seq_len(nchar(my_coding) / 3) %>% rep(each = 3) %>% rev)
  }
}


make_sequence_table <- function(my_gene) {
  cds = get_enst(my_gene) %>%
    get_cds_ranges

  coding = get_enst(my_gene) %>%
    get_coding_sequence %>%
    remove_stop_codon

  stopifnot(identical(width(cds) %>% sum,
                      nchar(coding)))

  data.table(chrom = get_chromosome(cds),
               pos = get_positions(cds),
               ref = orient_coding_sequence(cds, coding),
             aapos = orient_protein_pos(cds, coding))
}


make_bed <- function(my_table) {
  my_table[, .(chrom,
                          chromStart = pos - 1,
                          chromEnd = pos)]
}


write_bed <- function(my_bed, fname) {
  fwrite(my_bed,
         fname,
         sep = "\t",
         col.names = F)
}


run_liftover <- function(my_table) {
  fname_bed_hg38 = "bed_input.txt"
  fname_bed_hg19 = "bed_output.txt"
  fname_bed_unmapped = "bed_unmapped.txt"

  bed_hg38 = make_bed(my_table)
  write_bed(bed_hg38, fname_bed_hg38)

  cmd = paste(LIFTOVER_EXECUTABLE,
            fname_bed_hg38,
            LIFTOVER_CHAINFILE,
            fname_bed_hg19,
            fname_bed_unmapped,
            collapse = " ")

  # //// Runs in command line
  system(cmd)
  # /////////////////////////

  stopifnot(identical(file.size(fname_bed_unmapped), 0))

  bed_hg19 = fread(fname_bed_hg19)
  names(bed_hg19) = names(bed_hg38)

  my_table[, .(chrom,
               pos,
               pos_hg19 = bed_hg19$chromEnd,
               aapos,
               ref)]
}


# Run
# ===

# Read in data
idmapping = read_idmapping(UNIPROT_GENCODE_IDMAPPINGS)
gtf = read_gtf(GENCODE_GTF)
enst_coding = read_coding_translations(ENST_CODING)
uniprot_ids = read_uniprot_ids(UNIPROT_IDS)


# Select your gene here
g = "ACTA1"


# Build summary table, pass it through LiftOver
x = make_sequence_table(g)
x %<>% run_liftover


# Add columns for additional metadata
cols_idmapping = idmapping[enst == get_enst(g), .(hgnc = gene,
                                                        ensg,
                                                        transcript_id = enst,
                                                        assembly = "GRCh37")]

cols_uniprot_id = uniprot_ids[e == get_enst(g)]


# Reformat chr column, add genomic position column
x[, v_hg19 := paste0(chrom %>% gsub("chr", "", .),
                     ":",
                     pos_hg19,
                     ref)]


# Bind together
x = cbind(x, cols_idmapping, u)


# Rename and filter columns
x %<>% .[, .(uniprot_id = u,
             hgnc,
             ensg,
             transcript_id,
             protein_position = aapos,
             assembly,
             genomic_position = v_hg19,
             reference = ref)]


# Write to file
fwrite(x,
       paste0("positions_", g, ".txt"),
       sep = "\t")
