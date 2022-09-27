# LT 23/09/2022

library(tidyverse)

oe_file <- "gnomad/prototype/ENST00000366684_bp_oe.tsv"
align_file <- "structures/AF-P68133-F1-model_v2.align.tsv"

oe <- read_tsv(oe_file, col_types = "ccid")

# CDS length 1491 bp
print(nrow(oe))

# matches gnomADv2.1.1
print(sum(oe$exp))
print(sum(oe$obs))

hist(oe$exp)
plot(density(oe$exp), main = "expected missense per bp", xlab = "")

# note that loci with 0 expected missense are missing from oe
# add them

align <- read_tsv(align_file, col_type = "---ci-c-",
col_names = c("transcript", "codon", "locusR"))

align <- align %>% transmute(
        locus = sub('.$', '', locusR), codon)

print(tail(oe %>% right_join(align)))
exp_vect = oe %>% right_join(align) %>%
 transmute(exp = ifelse(is.na(exp), 0, exp)) %>%
 pull(exp)
plot(density(exp_vect))

# aggregate at codon level
oe_codon <- oe %>% right_join(align) %>% group_by(codon) %>%
summarise(
    obs = sum(obs, na.rm = TRUE),
    exp = sum(exp, na.rm = TRUE))

plot(density(oe_codon$exp), main = "expected missense per codon")
#hist(oe_codon$exp)

# save codon level oe file 
out_file <- "gnomad/prototype/ENST00000366684_oe.tsv"
write_tsv(oe_codon, out_file, col_names = FALSE)