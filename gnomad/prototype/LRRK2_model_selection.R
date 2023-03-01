# LT 22/12/2022
# experimental model selection on LRRK2

library(tidyverse)

GENE_ID <- "LRRK2"
TRANSECT_ID <- "ENST00000298910"
UNIPROT_ID <- "Q5S007"
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


getLL <- function(partition) {
    # partition is a list of sets of aa
    sum(sapply(partition, function(zone) {
        loe = oe %>% slice(zone)
        gamma = sum(loe$obs)/sum(loe$exp)
        nLL = -dpois(loe$obs, gamma * loe$exp, log = TRUE)
        sum(nLL)
    }))
}

#loop over all
nLLs = sapply(1:nrow(bubbles), function(i) {
    lsel = selected %>% add_row(slice(bubbles, i))
    lsel$bubble[[1]] <- setdiff(lsel$bubble[[1]], lsel$bubble[[2]])
    getLL(lsel$bubble)
})
print(min(nLLs))


# stop here for the sake of the demo: one round of model selection is 
# enough for LRRK2
# completing the round and starting the next one would involve
# get AIC and check if its lower than best AIC so far. if not stop
# update selected with new bubble
# remove bubble from pool candidate
# see example in ACTA1_model_selection.R

# likelihood ratio test
onezoneLL=getLL(selected$bubble)
twozoneLL = min(nLLs)
1-pchisq(2*(onezoneLL-twozoneLL), df=1)

# get the partition
bubbles[which.min(nLLs),]
ibest <- which.min(nLLs)

loe <- oe %>% slice(setdiff(1:n, bubbles$bubble[[ibest]]))
gamma <- sum(loe$obs)/sum(loe$exp)

#save result
res <- data.frame(index=1:n, gamma = 0)
res$gamma[setdiff(1:n, bubbles$bubble[[ibest]])] <- gamma
res$gamma[bubbles$bubble[[ibest]]] <- bubbles$oe[ibest]

gam <- round(res$gamma * 100)

# output a 2 columns file for iCn3D:
# first column: residue index
# second column:oe value on 0-100 range
write.table(gam,
    file = "gnomad/prototype/ENST00000298910_sel.tsv", sep = "\t",
    row.names = seq_len(length(gam)), col.names = FALSE, quote = FALSE
)
