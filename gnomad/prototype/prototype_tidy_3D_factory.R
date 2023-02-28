# LT 08/12/2022
# exploration of 3D factory with tidyverse

library(tidyverse)

GENE_ID = "ATP5F1B"
TRANSECT_ID = "ENST00000262030"
UNIPROT_ID = "P06576"
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

bubbles <- oe %>% pull(residue_index) %>% map(~dist_mat[,.])
    add_column(dist = dist_mat[20,])

    
     %>%
    arrange(dist) 

bubbles <- oe %>% mutate(dist = map(residue_index, ~ sort(dist_mat[,.])))

# prototype apply function to elements of list
test <- 
  tibble(
  x = c(1, 9, 5, 5),
  y = c(TRUE, FALSE, FALSE, FALSE),
  z = c("apple", "pear", "banana", "banana"),
  l = list(3:1, 4:1, 7:3, 7:3)
)
test %>% rowwise() %>% mutate(sorted = list(sort(l))) %>% pull(sorted)

# try aoe %>% add_column(dist = dist_mat[20,]) %>% arrange(dist) %>%
bubbles <- oe %>%
    mutate(zone = map_df(residue_index,
        ~ oe %>% add_column(dist = sort(dist_mat[,.])) %>%
            arrange(dist) %>%
            mutate(cum_exp = cumsum(exp), cum_obs = cumsum(obs), oe = cum_obs/cum_exp)))

mutate(cum_exp = cumsum(exp), cum_obs = cumsum(obs), oe = cum_obs/cum_exp) nested tible

bubbles <- oe %>%
    mutate(dist = map(residue_index, ~ sort(dist_mat[,.]))) %>%
    mutate(zone = list(oe)) %>%
    rowwise() %>% transmute(zone = list(zone %>% arrange(dist)))

bubbles <- oe %>%
    mutate(dist = map(residue_index, ~ dist_mat[,.])) %>%
    rowwise() %>%
    mutate(zone = list(
        oe %>% arrange(dist) %>%
        mutate(across(c(obs, exp), cumsum)) %>%
        mutate(oe_up = get_oe_ci_up(obs, exp)) %>%
        mutate(oe = obs / exp) %>%
        slice_min(oe_up)))
    
bubbles <- map_dfr(oe$residue_index,
    ~ arrange(oe, dist_mat[,.]) %>%
        mutate(
            index = 1:n,
            across(c(obs, exp), cumsum),
            oe = obs / exp,
            oe_up = get_oe_ci_up(obs, exp)) %>%
        slice_min(oe_up))
       
bubbles <- map_dfr(oe$residue_index, ~ {
    bubble = order(dist_mat[,.])
    arrange(oe, dist_mat[,.]) %>%
        mutate(
            index = 1:n,
            across(c(obs, exp), cumsum),
            oe = obs / exp,
            oe_up = get_oe_ci_up(obs, exp)) %>%
        slice_min(oe_up) %>%
        mutate(bubble = list(bubble[1:index])) %>%
        select(obs, exp, oe, oe_up, bubble)
    })

bubbles <- map_dfr(oe$residue_index, ~ {
    bubble = order(dist_mat[,.])
    slice(oe, bubble) %>%
        mutate(
            index = 1:n,
            across(c(obs, exp), cumsum),
            oe = obs / exp,
            oe_up = get_oe_ci_up(obs, exp)) %>%
        slice_min(oe_up) %>%
        mutate(bubble = list(sort(bubble[1:index]))) %>%
        rename(len = index) %>%
        select(obs, exp, oe, oe_up, len, bubble)
    }) %>% distinct


## prototype forward selection

# remove whole protein bubble if present
bubbles <- bubbles %>% filter(len < n)

# initialize tibble of selected bubbles whith whole protein bubble
selected <- oe %>% summarise(across(c(obs, exp), sum)) %>%
    mutate(
        oe = obs / exp,
        oe_up = get_oe_ci_up(obs, exp),
        len = n,
        bubble = list(1:n))

# get model corresponding to selected list
model <- selected %>%
    # build model with greedy adjudicator
    unnest(bubble) %>% rename(residue_index = bubble) %>%
    group_by(residue_index) %>% slice_min(oe) %>% select(residue_index, gamma = oe) %>% ungroup() %>%
    inner_join(oe) %>%
    # negative log-likelihood and AIC
    mutate(nLL = dpois(obs, gamma * exp))

best_AIC = 2 * sum(model$nLL) + 2 * 
        AIC = 2 * nLL + 2 * n_distinct(gamma)




# one round of forward selection
