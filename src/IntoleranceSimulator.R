# LT 4/05/2022


# Intolerance simulation module

# Use a statistical model to generate correlated random intolerance
# scores (o/e) on alphafold protein structures


## required packages
# for multivariate Gaussian RNG
require(MASS)

source('constraint_utils.R')

# protein structure file
# hard-coded for quick prototyping

prot_file <- "../structures/AF-P68133-F1-model_v2.pdb"

dist_mat = get_dist_mat(prot_file)

# calculate average distance between neighboring residues
aadist <- mean(diag(dist_mat[-nrow(dist_mat), -1]))

# spatial covariance function
# range:
l <- aadist * 5

# use a Gaussian kernel for prototyping
# sigma = exp(-0.5 * (dist_mat/l)^2)

# exponential kernel is less smooth
sigma <- exp(- dist_mat / l)


## check that sigma is positive definite
# tolerance needed because of numerical errors
tol <- 1e-6
if (min(eigen(sigma)$values) < -tol) stop("Spatial covariance matrix is not positive definite")

# generate a multivariate Gaussian
obs <- mvrnorm(1, mu = rep(0, nrow(sigma)), Sigma = sigma)


# use a copula to generate correlated multivariate beta observations

# hyper parameters: values matching the empirical distribution of o/e in gnomAD should be used
# here we fix an arbitrary beta to control the variance and calculate alpha such as the expected value is 0.5
beta <- 1
m <- 0.5
alpha <- beta * m / (1 - m)
obs <- qbeta(pnorm(obs), shape1 = alpha, shape2 = beta)

# output a 2 columns file for the 3D Factory:
# first column: residue index
# second column:simulated oe value
write.table(obs,
  file = "../structures/AF-P68133-F1-model_v2.simoe.tsv", sep = "\t",
  row.names = 1:nrow(sigma), col.names = FALSE, quote = FALSE
)

## Output for iCn3D (demo)
export_iCn3D(obs, "../structures/AF-P68133-F1-model_v2.simoe.iCn3D.tsv")