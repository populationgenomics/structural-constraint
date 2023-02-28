# LT 15/12/2022
# experimental model selection on ACTA1

library(tidyverse)

GENE_ID <- "ACTA1"
TRANSECT_ID <- "ENST00000366684"
UNIPROT_ID <- "P68133"
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


# tuning parameter
alpha <- 0.01

# tidyverse 3D Factory
bubbles <- map_dfr(oe$residue_index, ~ {
    bubble = order(dist_mat[,.])
    slice(oe, bubble) %>%
        mutate(
            index = 1:n,
            across(c(obs, exp), cumsum),
            oe = obs / exp,
            oe_up = get_oe_ci_up(obs, exp, alpha)) %>%
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
    mutate(nLL = dpois(obs, gamma * exp)) %>%
    # get AIC
    summarise(AIC = 2 * sum(nLL) + 2 * n_distinct(gamma))

#best_AIC = 2 * sum(model$nLL) + 2 * nrow(selected)


getLL <- function(partition) {
    # partition is a list of sets of aa
    sum(sapply(partition, function(zone) {
        loe = oe %>% slice(zone)
        gamma = sum(loe$obs)/sum(loe$exp)
        nLL = -dpois(loe$obs, gamma * loe$exp, log = TRUE)
        sum(nLL)
    }))
}

# choose a bubble to add 
#chosen <- bubbles %>% slice_min(oe) %>% slice_min(oe_up)

# should check its unique, I know it is here


# add to selected
#selected <- selected %>% add_row(chosen)

# remove duplicate residue from top zone
#selected$bubble[[1]] <- setdiff(selected$bubble[[1]], selected$bubble[[2]])

#loop over all
nLLs = sapply(1:nrow(bubbles), function(i) {
    lsel = selected %>% add_row(slice(bubbles, i))
    lsel$bubble[[1]] <- setdiff(lsel$bubble[[1]], lsel$bubble[[2]])
    getLL(lsel$bubble)
})
print(min(nLLs))

# likelihood ratio test
onezoneLL=getLL(selected$bubble)
twozoneLL = min(nLLs)
1-pchisq(2*(onezoneLL-twozoneLL), df=1)

# get the partition
bubbles[which.min(nLLs),]
ibest <- which.min(nLLs)

loe = oe %>% slice(setdiff(1:n, bubbles$bubble[[ibest]]))
        gamma = sum(loe$obs)/sum(loe$exp)
        nLL = -dpois(loe$obs, gamma * loe$exp, log = TRUE)

#save result
res = data.frame(index=1:n, gamma = 0)
res$gamma[setdiff(1:n, bubbles$bubble[[ibest]])] = gamma
res$gamma[bubbles$bubble[[ibest]]] = bubbles$oe[ibest]

gam <- round(res$gamma * 100)

# output a 2 columns file for iCn3D:
# first column: residue index
# second column:oe value on 0-100 range
write.table(gam,
    file = "gnomad/prototype/ENST00000366684_sel.tsv", sep = "\t",
    row.names = seq_len(length(oe)), col.names = FALSE, quote = FALSE
)



# all the rest bellow should be ignored 
# one round of forward selection

# do a loop for now
for (i in 1:nrow(bubbles)) {
    selected %>% add_row(slice(bubbles, i)) %>%
    unnest(bubble) %>% rename(residue_index = bubble) %>%
    group_by(residue_index) %>% slice_min(oe) %>% select(residue_index, gamma = oe) %>% ungroup() %>%
    inner_join(oe) %>%
    # negative log-likelihood and AIC
    mutate(nLL = dpois(obs, gamma * exp)) %>%
    # get AIC
    summarise(AIC = 2 * sum(nLL) + 2 * n_distinct(gamma))
}
bubbles %>% rowid_to_column("index") %>% rowwise() %>% mutate(
    model = list(
        selected %>% add_row(cur_data()) %>%
        unnest(bubble) %>% rename(residue_index = bubble) %>%
        group_by(residue_index) %>% slice_min(oe) %>% select(residue_index, gamma = oe) %>% ungroup() %>%
        inner_join(oe) %>%
        # negative log-likelihood and AIC
        mutate(nLL = dpois(obs, gamma * exp)) %>%
        # get AIC
        summarise(AIC = 2 * sum(nLL) + 2 * n_distinct(gamma))))

bubbles %>% rowid_to_column("index") %>% rowwise() %>% mutate(
    model = add_row(selected,cur_data())
        selected %>% add_row(cur_data())

# problems with cur_data(), try map
AICs = map(1:nrow(bubbles), ~ {
    i = .
    selected %>% add_row(slice(bubbles, i)) %>%
    unnest(bubble) %>% rename(residue_index = bubble) %>%
    group_by(residue_index) %>% slice_min(oe) %>% select(residue_index, gamma = oe) %>% ungroup() %>%
    inner_join(oe, by="residue_index") %>%
    # negative log-likelihood and AIC
    mutate(nLL = dpois(obs, gamma * exp)) %>%
    # get AIC
    summarise(AIC = 2 * sum(nLL) + 2 * n_distinct(gamma)) %>% pull(AIC)
})

AICs = sapply(1:nrow(bubbles), function(i) {
    selected %>% add_row(slice(bubbles, i)) %>%
    unnest(bubble) %>% rename(residue_index = bubble) %>%
    group_by(residue_index) %>% slice_min(oe) %>% select(residue_index, gamma = oe) %>% ungroup() %>%
    inner_join(oe, by="residue_index") %>%
    # negative log-likelihood and AIC
    mutate(nLL = dpois(obs, gamma * exp)) %>%
    # get AIC
    summarise(AIC = 2 * sum(nLL) + 2 * n_distinct(gamma)) %>% pull(AIC) 
})

bub <- bubbles %>% rowid_to_column("index")
 mutate(bub, AIC = map(index, ~ {
    add_row(selected, slice(bubbles,.)) %>%
    unnest(bubble) %>% rename(residue_index = bubble) %>%
    group_by(residue_index) %>% slice_min(oe) %>% select(residue_index, gamma = oe) %>% ungroup() %>%
    inner_join(oe, by="residue_index") %>%
    # negative log-likelihood and AIC
    mutate(nLL = dpois(obs, gamma * exp)) %>%
    # get AIC
    summarise(AIC = 2 * sum(nLL) + 2 * n_distinct(gamma)) %>% pull(AIC)
 }))

bub <- bubbles %>% rowid_to_column("index")
ms = mutate(bub, AIC = map(index, ~ {
    add_row(selected, slice(bubbles,.)) %>%
    unnest(bubble) %>% rename(residue_index = bubble) %>%
    group_by(residue_index) %>% slice_min(oe) %>% select(residue_index, gamma = oe) %>% ungroup() %>%
    inner_join(oe, by="residue_index)
}))
m