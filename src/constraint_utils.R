# LT 18/05/2022

# utility function for 3D missense constraint metric

## required packages
require(tidyverse)
# for .PDB structures
require(bio3d)

# set seed for reproducibility
set.seed(9011971)

## pairwise distances between residues based on their carbon alpha
get_dist_mat <- function(pdb_file) {
  pdb <- read.pdb(pdb_file)

  # calculate distance between residues based on their carbon alpha
  dist_mat <- dm(pdb, inds = "calpha")

  # dist_mat is upper triangular with empty diagonal:
  # add diagonal, fill lower triangle, convert to matrix
  dist_mat <- as.matrix(dist_mat)
  dist_mat[is.na(dist_mat)] <- 0
  dist_mat <- dist_mat + t(dist_mat)

  dist_mat
}

## export a dataset containing intolerance scores in a format compatible with
## iCn3D

export_iCn3D <- function(oe_data, file) {
  if ("numeric" %in% class(oe_data)) {
    oe <- oe_data
  } else {
    oe <- oe_data[, "oe"]
  }

  # set score on 0-100 range (needed by iCn3D)
  oe <- round(oe * 100)

  # output a 2 columns file for iCn3D:
  # first column: residue index
  # second column:simulated oe value on 0-100 range
  write.table(oe,
    file = file, sep = "\t",
    row.names = seq_len(length(oe)), col.names = FALSE, quote = FALSE
  )
}

## calculate upper bound of the o/e confidence interval
get_oe_ci_up <- function(observed, expected, alpha = 0.05) {
  return(qchisq((1 - alpha / 2), 2 * (observed + 1)) / 2 / expected)
}
