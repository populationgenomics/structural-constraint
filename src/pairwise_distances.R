# pairwise_distances.R
# M. Silk

# Calculate pairwise distances between
# residues within the AF2 structures


suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(magrittr))
suppressPackageStartupMessages(library(purrr))
suppressPackageStartupMessages(library(bio3d))


args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 2)


file_in <- args[1]
file_out <- args[2]
stopifnot(file.exists(file_in))


read_structure <- function(filename) {
  read.pdb(filename) %>%
    return
}


extract_atom_table <- function(pdb) {
  pdb[["atom"]] %>%
    as.data.table %>%
    return
}


calc_atom_dists <- function(atoms, method = "euclidean") {
  atoms[, .(x, y, z)] %>%
    dist(method = method) %>%
    as.matrix %>%
    round(3) %>%
    as.data.table %>%
    return
}


# Summarise atom distances by assigning residues
# to the atoms, taking minimum distance between
# each RESIDUE to every atom, flip table and
# repeat, taking min dist from each RESIDUE to
# RESIDUE.
min_resi_dists <- function(atom_dists, atoms) {
  cbind(atoms[, .(resno)], atom_dists) %>%
    .[, map(.SD, min), resno] %>%
    .[, resno := NULL] %>%
    t %>%
    as.data.table %>%
    cbind(atoms[, .(resno)], .) %>%
    .[, map(.SD, min), resno] %>%
    .[, resno := NULL] %>%
    return
}


# Bind resid, chain, resno columns to the dist table
make_dist_table <- function(resi_dists, atoms) {
  cbind(atoms[, .(resid, chain, resno)] %>% .[!duplicated(.)],
        resi_dists) %>%
    return
}


write_output_table <- function(dist_table, file) {
  fwrite(dist_table,
         file,
         sep = "\t")
}


# Run
# atoms table is used in several steps and so is read
# into memory as an object, and also passed to most of the
# subsequent functions.
atoms <- read_structure(file_in) %>%
  extract_atom_table

calc_atom_dists(atoms) %>%
  min_resi_dists(atoms) %>%
  make_dist_table(atoms) %>%
  write_output_table(file_out)
