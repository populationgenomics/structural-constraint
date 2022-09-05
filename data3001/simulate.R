# LT 5/09/2022


# Simulate data for DATA3001

# Use a statistical model to generate correlated random intolerance
# scores (o/e) on alphafold protein structures


## source utils
source("src/constraint_utils.R")

# protein structure files for the 10 genes initial panel
pdb_files <- c(
  "AF-O43526-F1-model_v3.pdb",
  "AF-P01009-F1-model_v3.pdb",
  "AF-P06576-F1-model_v3.pdb",
  "AF-P17181-F1-model_v3.pdb",
  "AF-P60484-F1-model_v3.pdb",
  "AF-P68133-F1-model_v3.pdb",
  "AF-P69905-F1-model_v3.pdb",
  "AF-Q5S007-F1-model_v3.pdb",
  "AF-Q86VV8-F1-model_v3.pdb",
  "AF-Q969H0-F1-model_v3.pdb")

for (pdb_file in pdb_files) {

  # model name
  model_name <- tools::file_path_sans_ext(pdb_file)

  pdb_file <- paste("structures/", pdb_file, sep = "")

  # read PDB file
  pdb <- read.pdb(pdb_file)

  ## Select calpha atoms
  calpha <- atom.select(pdb, "calpha")

  ## Trim XYZ to get coordinates of carbon alphas
  my_xyz <- trim(pdb$xyz, col.inds = calpha$xyz)

  # convert to coordinates matrix
  xyz_mat <- matrix(as.vector(my_xyz), ncol = 3, byrow = TRUE)

  # number of amino acids in protein
  n <- nrow(xyz_mat)

  # create data structure for this protein
  sim_data <- data.frame(index = seq_len(n),
    x = xyz_mat[, 1], y = xyz_mat[, 2], z = xyz_mat[, 3])

  # simulate gamma
  sim_data$gamma <- simulate_gamma(pdb_file)

  # simulate expected number of missense per residue
  # use gnomAD v2.2.1 parameters
  sim_data$exp <- rnorm(n, mean = 0.558, sd = 0.0519)

  # simulate observed number of missense per residue
  sim_data$obs <- rpois(n, sim_data$exp * sim_data$gamma)

  # write output
  write.table(sim_data,
              file = paste("data3001/", model_name, ".tsv", sep = ""),
              sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE
  )
}
