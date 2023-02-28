# LT 20/10/2022

library(tidyverse)

GENE_ID <- "ATP5F1B"
TRANSECT_ID <- "ENST00000262030"
UNIPROT_ID <- "P06576"

GENE_ID <- "LRRK2"
TRANSECT_ID <- "ENST00000298910"
UNIPROT_ID <- "Q5S007"

oe_file <- sprintf("gnomad/prototype/%s_bp_oe.tsv", TRANSECT_ID)
align_file <- sprintf("structures/positions_%s.txt", GENE_ID)

oe <- read_tsv(oe_file, col_types = "ccdd")

# CDS length
print(nrow(oe))

# check if matches gnomADv2.1.1
print(sum(oe$exp))
print(sum(oe$obs))

hist(oe$exp)
plot(density(oe$exp), main = "expected missense per bp", xlab = "")

# note that loci with 0 expected missense are missing from oe
# add them

align <- read_tsv(align_file, col_type = "---ci-c-",
col_names = TRUE)

align <- align %>% transmute(
        locus = genomic_position, codon =protein_position)

print(tail(oe %>% right_join(align)))
exp_vect <- oe %>% right_join(align) %>%
 transmute(exp = ifelse(is.na(exp), 0, exp)) %>%
 pull(exp)
plot(density(exp_vect), main = "distribution of expected missense per locus",
xlab = "number of expected missense", cex = 1.2)

# aggregate at codon level
oe_codon <- oe %>% right_join(align) %>% group_by(codon) %>%
summarise(
    obs = sum(obs, na.rm = TRUE),
    exp = sum(exp, na.rm = TRUE))

plot(density(oe_codon$exp), main = "distribution of expected missense per codon",
xlab = "number of expected missense")
#hist(oe_codon$exp)

# save codon level oe file 
out_file <- sprintf("gnomad/prototype/%s_oe.tsv", TRANSECT_ID)
write_tsv(oe_codon, out_file, col_names = FALSE)

source("src/constraint_utils.R")

# protein structure file. Alphafold v3 model
prot_file <- sprintf("structures/AF-%s-F1-model_v3.pdb", UNIPROT_ID)

# o/e file, out of the GenomicSimulator or a genomic database
oe_file <- sprintf("gnomad/prototype/%s_oe.tsv", TRANSECT_ID)

oe <- read_tsv(oe_file, col_names = c("residue_index", "obs", "exp"))

# calculate pairwise distances between residues based on their carbon alpha
dist_mat <- get_dist_mat(prot_file)

# number of residues in protein
n <- nrow(dist_mat)


# plot sphere scores

oe2 <- oe %>% add_column(dist = dist_mat[20,]) %>% arrange(dist) %>%
mutate(cum_exp = cumsum(exp), cum_obs = cumsum(obs), oe = cum_obs/cum_exp)

alphas <- c(0.4, 0.3, 0.2, 0.1, 0.05, 0.01, 0.001, 0.0001, 0.00001)
oe_up <- sapply(alphas, function(alpha) get_oe_ci_up(oe2$cum_obs, oe2$cum_exp, alpha))
colnames(oe_up) <- alphas

oe2 <- oe2 %>% bind_cols(oe_up)

fig = ggplot(oe2) + geom_point(aes(x = dist, y = oe), size = 2) + theme_classic(base_size=15) +
geom_line(aes(x = dist, y = `0.2`), lwd = 2, col = 2) + 
geom_line(aes(x = dist, y = `0.05`), lwd = 2, col = 3) +
geom_line(aes(x = dist, y = `0.001`), lwd = 2, col = 4) +
scale_y_log10()

plot(fig)

## calculate intolerance metric of a region
get_region_metric <- function(i_region, alpha, regions, oe_data) {
  # indices of the residues belonging to this region
  this_region <- regions[seq_len(i_region)]

  # total number of missense variants observed in the region
  total_obs <- sum(oe_data$obs[this_region])

  # total number of missense variants expected in the region
  total_exp <- sum(oe_data$exp[this_region])

  # upper bound of o/e
  oe_upper <- get_oe_ci_up(total_obs, total_exp, alpha)

  c(obs = total_obs, exp = total_exp, oe = total_obs / total_exp, oe_upper = oe_upper)
}

# get 3D region and its metrics for each residue
regions_list <- lapply(seq_len(n), function(i_residue) {

  # get the corresponding row of the distance matrix
  row <- dist_mat[i_residue, ]

  # grow the sphere = order by increasing distance
  this_aa_regions <- seq_len(n)[order(row)]

  # calculate the upper bound of the o/e confidence interval for each region
  regions_metric <- sapply(seq_len(n), function(i) {
    get_region_metric(i, 0.0001, this_aa_regions, oe_data = oe)
  })

  # find index of best region: that with the lowest upper bound of o/e CI
  i_region <- which.min(regions_metric["oe_upper", ])

  # vector of the indices of the residues belonging to the region
  region <- this_aa_regions[seq_len(i_region)]

  # return data structure
  list(metric = regions_metric[, i_region], region = region)
})

# Adjudicator

# creates a partition of the protein into regions of intolerance
# assigns to each residue a single intolerance score
# prototype: greedy algorithm: discussion with Sabrina and Leo, 28/04/2022

# sort list of regions by increasing intolerance score and create map of residue to region

# intolerance scores
oes <- sapply(regions_list, function(one_region) one_region$metric["oe"])

# sort list
regions_list <- regions_list[order(oes)]

# create empty map from residue to region
residue_to_region <- cbind(residue_index = seq_len(n), region_index = NA)

# fill the map, going down the list of ordered regions
for (i_region in seq_len(n)) {
  this_region <- regions_list[[i_region]]
  # residues in this region that do not have an intolerance score yet
  to_score <- which(is.na(residue_to_region[this_region$region, "region_index"]))

  # assign the regions's intolerance score to these residues
  residue_to_region[this_region$region, "region_index"][to_score] <- i_region

  # end if all residues have been scored
  if (all(!is.na(residue_to_region[, "region_index"]))) break
}

## decorate map with intolerance metric

intolerance_scores <- t(apply(residue_to_region, 1, function(one_row) {
  # get metric information from regions_list
  metric <- regions_list[[one_row["region_index"]]]$metric
  c(one_row, metric)
}))

#
# save results

# intolerance score and region per residue
write.table(intolerance_scores,
  file = sprintf("gnomad/prototype/%s_intol.tsv", TRANSECT_ID), sep = "\t",
  row.names = FALSE, col.names = TRUE
)

# for visualisation in iCn3D
intolerance_scores[,"oe"] = intolerance_scores[, "oe"]/max(intolerance_scores[,"oe"])
export_iCn3D(intolerance_scores, sprintf("gnomad/prototype/%s_intol_iCn3D.tsv", TRANSECT_ID))
