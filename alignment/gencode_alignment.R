# gencode_alignment.R
# M. Silk

# Align the GENCODE transcripts to the AlphaFold2
# protein tertiary structures

# If running on gnomAD v2.1.1, pipeline is set up
# to use GENCODE v19 / GRCh37 coordinates

library(data.table)
library(magrittr)
library(purrr)
library(bio3d)


# Define directory for AF2 structures
dir_af2 <- "~/data/af2_human_v3/"


# Read in the UniProt idmapping table
idmapping_names <- c(
  "UniProtKB-AC",
  "UniProtKB-ID",
  "GeneID (EntrezGene)",
  "RefSeq",
  "GI",
  "PDB",
  "GO",
  "UniRef100",
  "UniRef90",
  "UniRef50",
  "UniParc",
  "PIR",
  "NCBI-taxon",
  "MIM",
  "UniGene",
  "PubMed",
  "EMBL",
  "EMBL-CDS",
  "Ensembl",
  "Ensembl_TRS",
  "Ensembl_PRO",
  "Additional PubMed"
)
idmapping <- fread("alignment/HUMAN_9606_idmapping_selected.tab.gz", col.names = idmapping_names)


# Read in filenames of AF2 structures and their UniProt IDs
# Remove any AF2's with multiple fragments
# Note that with current release, there is only one structure
# per UniProt
af2_files <- list.files(path = dir_af2, pattern = "*.pdb.gz")
af2 = data.table(af2_file = af2_files)
af2[, uniprot := af2_file %>% strsplit("-") %>% map(2) %>% unlist]
uniprots_with_multifrags <- af2[grep("-F2-", af2_file, fixed = TRUE), uniprot] %>%
  unique
af2 <- af2[!uniprot %in% uniprots_with_multifrags]


# Filter UniProt accessions to those with an AF2 structure
idmapping_af2 <- idmapping[`UniProtKB-AC` %in% af2$uniprot]


# Simplify table to a list of UniProts and their listed ENSGs
idmapping_af2_ensgs <- idmapping_af2[
  !is.na(Ensembl), .(
    uniprot = `UniProtKB-AC`,
    ensg = strsplit(Ensembl, "; ", fixed = TRUE)[[1]]
  ),
  seq_len(nrow(idmapping_af2))
] %>%
  .[, seq_len := NULL]


# Read in the GENCODE v19 translations and remove version numbers to
# ENSTs and ENSGs
gencode_translations <- fread("alignment/gencode_translations.txt.gz")
gencode_translations[, ensg := strsplit(ensg, ".", fixed = TRUE) %>% map(1) %>% unlist]
gencode_translations[, enst := strsplit(enst, ".", fixed = TRUE) %>% map(1) %>% unlist]


# Bind uniprot ID's to gencode_translations
setkey(gencode_translations, ensg)
setkey(idmapping_af2_ensgs, ensg)
gencode_translations_uniprots = idmapping_af2_ensgs[gencode_translations]


# Read a PDB sequence from UniProt ID
read_structure_peptide <- function(my_uniprot) {
  fname = paste0(dir_af2, af2[uniprot == my_uniprot]$af2_file)
  read.pdb(fname) %>%
    pdbseq() %>%
    paste(collapse = "")
}


# Match GENCODE sequence to a protein sequence
check_match <- function(gencode_seq, af2_seq) {
  ifelse(gencode_seq == af2_seq, TRUE, FALSE)
}


# Match sequence to AF2 structure
gencode_translations_uniprots[!is.na(uniprot),
                              match := map2(sequence, read_structure_peptide(uniprot), check_match),
                              by = seq_len(nrow(gencode_translations_uniprots[!is.na(uniprot)]))]


# Write to file
fwrite(gencode_translations_uniprots, "alignment/gencode_translations_matched.txt.gz", sep = "\t")
