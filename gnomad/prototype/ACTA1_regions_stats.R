# LT 27/09/2022

# analyse prototype results for ACTA1

source("src/constraint_utils.R")

# protein structure file
prot_file <- "structures/AF-P68133-F1-model_v2.pdb"

# o/e file, out of the GenomicSimulator or a genomic database
oe_file <- "gnomad/prototype/ENST00000366684_oe.tsv"

oe <- read_tsv(oe_file, col_names = c("residue_index", "obs", "exp"))

# calculate pairwise distances between residues based on their carbon alpha
dist_mat <- get_dist_mat(prot_file)

# number of residues in protein
n <- nrow(dist_mat)

##################
# 3D Factory

oe2 <- oe %>% add_column(dist = dist_mat[2,]) %>% arrange(dist) %>%
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
# plot region

  # fig = ggplot(dat) + geom_line(aes(x = octave, y = pred),lwd=1) +
  #   geom_point(aes(x = octave, y = obs), size=2, alpha=1) +
  #   theme_classic(base_size=15) + theme(legend.position = 'top') + scale_y_sqrt() +
  #   ggtitle(dat$gene) +
  #   ylab('Proportion of variants')

# get 3D region and its metrics for each residue
regions_list <- lapply(seq_len(n), function(i_residue) {

  # build a spheres object
  sphere <- list(center = i_residue)

  # get the corresponding row of the distance matrix
  row <- dist_mat[i_residue, ]

  # grow the sphere = order by increasing distance
  sphere$indices <- order(row)


  # order oe data by distance
  oe_data <- oe[sphere$indices,]
  # calculate the upper bound of the o/e confidence interval for each region
  regions_metric <- sapply(seq_len(n), function(i) {
    get_region_metric(i, 0.0001, this_aa_spheres, oe_data = oe)
  })

  # find index of best region: that with the lowest upper bound of o/e CI
  i_region <- which.min(regions_metric["oe_upper", ])

  # vector of the indices of the residues belonging to the region
  region <- this_aa_spheres[seq_len(i_region)]

  # return data structure
  list(center = i_residue, dist = row[this_aa_spheres], spheres = this_aa_spheres, metric = regions_metric, region_end = i_region)
})

################
# Adjudicator

# creates a partition of the protein into regions of intolerance
# assigns to each residue a single intolerance score
# prototype: greedy algorithm: discussion with Sabrina and Leo, 28/04/2022

## sort list of regions by increasing intolerance score and create map of residue to region

# # intolerance scores
# oes <- sapply(regions_list, function(one_region) one_region$metric["oe"])

# # sort list
# regions_list <- regions_list[order(oes)]

# # create empty map from residue to region
# residue_to_region <- cbind(residue_index = seq_len(n), region_index = NA)

# # fill the map, going down the list of ordered regions
# for (i_region in seq_len(n)) {
#   this_region <- regions_list[[i_region]]

#   # residues in this region that do not have an intolerance score yet
#   to_score <- which(is.na(residue_to_region[this_region$region, "region_index"]))

#   # assign the regions's intolerance score to these residues
#   residue_to_region[this_region$region, "region_index"][to_score] <- i_region

#   # end if all residues have been scored
#   if (all(!is.na(residue_to_region[, "region_index"]))) break
# }

# ## decorate map with intolerance metric

# intolerance_scores <- t(apply(residue_to_region, 1, function(one_row) {
#   # get metric information from regions_list
#   metric <- regions_list[[one_row["region_index"]]]$metric
#   c(one_row, metric)
# }))

#
## save results

# intolerance score and region per residue
# write.table(intolerance_scores,
#   file = "gnomad/prototype/ENST00000366684_intol.tsv", sep = "\t",
#   row.names = FALSE, col.names = TRUE
# )

# # for visualisation in iCn3D
# intolerance_scores[,"oe"] = intolerance_scores[, "oe"]/max(intolerance_scores[,"oe"])
# export_iCn3D(intolerance_scores, "gnomad/prototype/ENST00000366684_intol_iCn3D.tsv")
