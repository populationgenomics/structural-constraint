# LT 18/05/2022

# utility function for 3D missense constraint metric

## required packages
# for .PDB structures
require(bio3d)
# for multivariate Gaussian RNG
require(MASS)

# set seed for reproducibility
set.seed(9011971)

## pairwise distances between residues based on their carbon alpha
get_dist_mat <- function(pdb_file) {
  pdb <- read.pdb(pdb_file)

  # calculate distance between residues based on their carbon alpha
  dist_mat <- dm(pdb, inds = "calpha")

  # dist_mat is upper triangular with empty diagonal:
  # convert to matrix, add diagonal, fill lower triangle
  dist_mat <- as.matrix(dist_mat)
  dist_mat[is.na(dist_mat)] <- 0
  dist_mat <- dist_mat + t(dist_mat)

  dist_mat
}

## export a dataset containing intolerance scores in a format compatible with
## iCn3D

export_iCn3D <- function(oe_data, file) { # nolint
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
    row.names = seq_len(nrow(sigma)), col.names = FALSE, quote = FALSE
  )
}

# simulate gamma, the true unobserved proportional depletion in missense mutation # nolint
# using a spatial correlation model
# and hyperparameters similar to observed distributions in gnomAD v2.1
simulate_gamma <- function(prot_file) {

  dist_mat <- get_dist_mat(prot_file)

  # calculate average distance between neighboring residues
  aadist <- mean(diag(dist_mat[-nrow(dist_mat), -1]))

  # spatial covariance function
  # range:
  l <- aadist * 5

  # use an exponential kernel for prototyping, it is less smooth than a Guassian
  sigma <- exp(-dist_mat / l)

  ## check that sigma is positive definite
  # tolerance needed because of numerical errors
  tol <- 1e-6
  if (min(eigen(sigma)$values) < -tol)
  stop("Spatial covariance matrix is not positive definite")

  # generate a multivariate Gaussian
  obs <- mvrnorm(1, mu = rep(0, nrow(sigma)), Sigma = sigma)

  # use a copula to generate correlated multivariate beta observations

  # hyper parameters:
  # values matching the empirical distribution of o/e in gnomAD should be used
  # here we fix an arbitrary beta to control the variance and calculate alpha
  # such that the expected value is 0.5
  beta <- 1
  m <- 0.5
  alpha <- beta * m / (1 - m)
  gamma <- qbeta(pnorm(obs), shape1 = alpha, shape2 = beta)

  return(gamma)
}