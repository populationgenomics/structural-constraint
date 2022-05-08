# LT 8/05/2022


# Power calculator

# a collection of code bits to calculate statistics of intolerance
# in gnomAD v2.1 and gnomAD v4


### estimate the distribution of the number of expected missense mutations per residue
## using the missense data at gene level from gnomAD v2.1

# gnomAD v2.1 intolerance data at gene level
cons = read.delim('../structures/gnomad.v2.1.1.lof_metrics.by_gene.txt.bgz')
exp_aa = cons$exp_mis / cons$cds_length * 3

summary(exp_aa)
hist(exp_aa)
# long tail to the left

# look at 10 genes with very low exp_aa
cons[order(exp_aa)[1:10],]

nrow(cons)
boxplot(exp_aa)

# remove outliers
outliers = boxplot(exp_aa, plot=FALSE)$out
clean_exp_aa = exp_aa[-which(exp_aa %in% outliers)]

# replot
hist(clean_exp_aa)
# normal approximation is OK here

# this gives us rough estimates to simulate expected missense per aa
expaa_mean = mean(clean_exp_aa)
expaa_sd = sd(clean_exp_aa)


## a rough power analysis to come below
