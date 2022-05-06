# LT 4/05/2022


# Intolerance simulation module

# Use a statistical model to generate correlated random intolerance
# scores (o/e) on alphafold protein structures


## required packages
require(bio3d)

# for multivariate Gaussian RNG
require(MASS)

# protein structure file

# hard-coded for quick prototyping

protfile = 'structures/AF-P68133-F1-model_v2.pdb'

pdb = read.pdb(protfile)

# calculate residues pairwise distance matrix
distmat = dm(pdb)

# distmat is upper triangular with empty diagonal: add diagonal, fill lower triangle, convert to matrix
dismat = as.matrix(distmat)
distmat[is.na(distmat)] = 0
distmat = distmat + t(distmat)

# spatial covariance function: use a Gaussian kernel for prototyping

l=1
sigma = exp(-0.5 * (distmat/l)^2)

# exponential kernel
l=1
sigma = exp(-0.5 * distmat / l)


# check that sigma is positive definite
if(min(eigen(sigma)$values) < 0) stop('Spatial covariance matrix is not positive definite')

# generate a multivariate Gaussian 
obs = rmvnorm(1, Sigma=sigma)


# use a copula to generate correlated multivariate beta observations

# hyper parameters: values matching the empirical distribution of o/e in gnomAD should be used
# here we fix an arbitrary beta to control the variance and calculate alpha such as the expected value is 0.4
beta = 5
m = 0.4
alpha = beta * m / (1-m)
obs = qbeta(pnorm(obs), shape1 = alpha, shape2 = beta)


# set score on 0-100 range
obs = round(obs * 100)

# output a 2 columns file for iCn3D
write.table(t(obs), file='structures/AF-P68133-F1-model_v2.simoe.tsv', sep='\t',
            row.names=1:nrow(sigma), col.names=FALSE, quote=FALSE)

