# LT 17/05/2022

# prototype: 3D Factory AND Adjudicator are in the same source file to avoid
# implementation of the interface between them

# 3D Factory component

# the 3D Factory generates a list of 3D regions and calculates their intolerance score

# First prototype
# creates 3D regions of variable length using upper bound of CI
# each region is built by growing a sphere centered on the carbon-alpha of
# a specific residue until the whole protein is included
# a criteria (based on the upper bound) is used to decide the length of the region

# file names and parameters are hard-coded for this prototype

source('constraint_utils.R')

## calculate upper bound of the o/e confidence interval
get_oe_ci_up = function (observed, expected, alpha = 0.05) {
  
  return(qchisq((1-alpha/2), 2*(observed + 1))/2/expected)
  
}

## calculate intolerance metric of a region
get_region_metric = function (i_region, alpha, regions, oe_data)
{
  # indices of the residues belonging to this region
  this_region = regions[seq_len(i_region)]
  
  # total number of missense variants observed in the region
  total_obs = sum(oe_data$obs[this_region])
  
  # total number of missense variants expected in the region
  total_exp = sum(oe_data$exp[this_region]) 
  
  # upper bound of o/e
  oe_upper = get_oe_ci_up(total_obs, total_exp, alpha)
  
  c(obs = total_obs, exp = total_exp, oe = total_obs/total_exp, oe_upper = oe_upper)
}

### main
## inputs 
# protein structure file
prot_file <- "../structures/AF-P68133-F1-model_v2.pdb"

# o/e file, out of the GenomicSimulator or a genomic database
oe_file = "../structures/AF-P68133-F1-model_v2.simmut.tsv"

oe <- read.table(file = oe_file, sep = "\t",
  col.names = c("residue_index", "obs", 'exp')
)

# calculate pairwise distances between residues based on their carbon alpha
dist_mat <- get_dist_mat(prot_file)

# number of residues in protein
n =nrow(dist_mat)

##################
# 3D Factory

# get 3D region and its metrics for each residue
regions_list = lapply(seq_len(n), function(i_residue) {
  
  # get the corresponding row of the distance matrix
  row = dist_mat[i_residue,]
  
  # grow the sphere = order by increasing distance
  this_aa_regions = seq_len(n)[order(row)]
  
  # calculate the upper bound of the o/e confidence interval for each region
  regions_metric = sapply(seq_len(n), function(i) { get_region_metric(i, 0.2, this_aa_regions, oe_data = oe) })
  
  # find index of best region: that with the lowest upper bound of o/e CI
  i_region = which.min(regions_metric['oe_upper',])
  
  # vector of the indices of the residues belonging to the region
  region = this_aa_regions[seq_len(i_region)]
  
  # return data structure
  list( metric = regions_metric[, i_region], region = region)
})

################
# Adjudicator

# creates a partition of the protein into regions of intolerance
# assigns to each residue a single intolerance score
# prototype: greedy algorithm: discussion with Sabrina and Leo, 28/04/2022

## sort list of regions by increasing intolerance score and create map of residue to region

# intolerance scores
oes = sapply(regions_list, function(one_region) one_region$metric['oe'])

# sort list
regions_list = regions_list[order(oes)]

# create empty map from residue to region
residue_to_region = cbind(residue_index = seq_len(n), region_index = NA)

# fill the map, going down the list of ordered regions
for (i_region in seq_len(n)) {
  this_region = regions_list[[i_region]]
  
  # residues in this region that do not have an intolerance score yet
  to_score = which(is.na(residue_to_region[this_region$region, 'region_index']))
  
  # assign the regions's intolerance score to these residues
  residue_to_region[this_region$region, 'region_index'][to_score] = i_region
  
  # end if all residues have been scored
  if (all(!is.na(residue_to_region[,'region_index']))) break
}

## decorate map with intolerance metric

intolerance_scores = t(apply(residue_to_region, 1, function(one_row) {
  # get metric information from regions_list
  metric = regions_list[[one_row['region_index']]]$metric
  c(one_row, metric) }))
  
#
## save results

# intolerance score and region per residue
write.table(intolerance_scores,
            file = "../structures/AF-P68133-F1-model_v2.intol.tsv", sep = "\t",
            row.names = FALSE, col.names = TRUE)

# for visualisation in iCn3D
export_iCn3D(intolerance_scores, "../structures/AF-P68133-F1-model_v2.intol.iCn3D.tsv")