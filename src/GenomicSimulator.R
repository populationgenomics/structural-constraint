# LT 16/05/2022


# Genomic Simulator component

# the genomic simulator component generates simulated data as if extracted from a large genomic database:
# counts of observed missense variants per residue
# expected number of missense variants

# First prototype
# simulate expected from a Gaussian distribution with parameters matching roughly gnomAD v2.2.1 (see PowerCalculator)
# Simulate counts from Poisson distribution with rate = expected
# all the file names and parameters are hard-coded for this prototype


# set seed for reproducibility
set.seed(9011971)

# open file with simulated o/e scores
oeTruth <- read.table(
  file = "../structures/AF-P68133-F1-model_v2.simoe.tsv", sep = "\t",
  col.names = c("residue_index", "oe")
)

# number of residues in protein
n <- nrow(oeTruth)

# simulate expected number of missense per residue
# use gnomAD v2.2.1 parameters

expected <- rnorm(n, mean = 0.558, sd = 0.0519)

# simulate observed number of missense per residue
observed <- rpois(n, expected * oeTruth$oe)

sim <- data.frame(oeTruth$residue_index, obs = observed, exp = expected)

write.table(sim,
  file = "../structures/AF-P68133-F1-model_v2.simmut.tsv", sep = "\t",
  row.names = FALSE, col.names = FALSE, quote = FALSE
)
