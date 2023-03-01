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
    bubble <- order(dist_mat[, .])
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

# functions for model selection
# calculate negative log-likelihood of a partition 
getLL <- function(partition) {
    # partition is a list of sets of aa covering the whole protein
    sum(sapply(partition, function(zone) {
        loe <- oe %>% slice(zone)
        gamma <- sum(loe$obs) / sum(loe$exp)
        nLL <- -dpois(loe$obs, gamma * loe$exp, log = TRUE)
        sum(nLL)
    }))
}

# calculate AIC of a partition
getAIC <- function(partition) {
    2 * sum(sapply(partition, length) > 0) + 2 * getLL(partition)
}

# remove whole protein bubble if present
bubbles <- bubbles %>% filter(len < n)

# initialize tibble of selected bubbles whith whole protein bubble
# this special bubble plays the role of "all the rest"
selected <- oe %>% summarise(across(c(obs, exp), sum)) %>%
    mutate(
        oe = obs / exp,
        oe_up = get_oe_ci_up(obs, exp),
        len = n,
        bubble = list(1:n))

# null model: only one zone ocvering the whole protein
null_model <- selected

# initialize best AIC
best_AIC <- getAIC(selected$bubble)

repeat {
    # perform one round of model selection

    cat("Starting new round of model selection\n")
    cat("best model df:", nrow(selected), " nLL:", getLL(selected$bubble), " AIC:", best_AIC, "\n")

    # add each bubble one by one and evaluate likelihood
    nLLs <- sapply(1:nrow(bubbles), function(i) {
        lsel <- selected %>% add_row(slice(bubbles, i))
        # remove what we added from the "all the rest" bubble
        lsel$bubble[[1]] <- setdiff(lsel$bubble[[1]], lsel$bubble[[nrow(lsel)]])
        getLL(lsel$bubble)
    })

    ## get AIC of best candidate model
    # index of best bubble this round
    ibest <- which.min(nLLs)[1]
    best_bubble <- slice(bubbles, ibest)
    candidate_model <- selected %>% add_row(slice(bubbles, ibest))
    # remove what we added from the "all the rest" bubble
    candidate_model$bubble[[1]] <- setdiff(candidate_model$bubble[[1]], candidate_model$bubble[[nrow(candidate_model)]])
    candidate_AIC <- getAIC(candidate_model$bubble)

    # stop if candidate model is not better
    if (candidate_AIC >= best_AIC) break

    # update best model
    best_AIC <- candidate_AIC
    selected <- candidate_model

    # update bubble list
    bubbles <- bubbles %>%
        mutate(
            bubble = map(bubble, function(zone) setdiff(zone, best_bubble$bubble[[1]])),
            len = map_int(bubble, length)) %>%
        filter(len > 0)
    if (nrow(bubbles) == 0) break

}

# likelihood ratio test comparing best model with null model
onezoneLL <- getLL(null_model$bubble)
bestLL <- getLL(selected$bubble)
p_val <- 1 - pchisq(2 * (onezoneLL - bestLL), df = nrow(selected) - 1)


#save result
res <- selected %>% select(bubble) %>%
    mutate(gamma = map_dbl(bubble, function(zone) {
        loe <- oe %>% slice(zone)
        gamma <- sum(loe$obs) / sum(loe$exp)
    })) %>%
    unnest(bubble) %>% rename(residue_index = bubble) %>%
    arrange(residue_index)

gam <- round(res$gamma * 100)

# output a 2 columns file for iCn3D:
# first column: residue index
# second column:oe value on 0-100 range
write.table(gam,
    file = "gnomad/prototype/ENST00000366684_sel.tsv", sep = "\t",
    row.names = seq_len(length(oe)), col.names = FALSE, quote = FALSE
)
