library(qtl2)
library(optparse)
library(dplyr)

option_list <- list(
  make_option(c("-j", "--json"), type = "character", default = NULL),
  make_option(c("-p", "--phenotype_file"), type = "character", default = NULL),
  make_option(c("-o", "--output_file"), type = "character", default = NULL),
  make_option(c("-m", "--mutation_type", type = "character", default = NULL)))

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

mutation_type <- gsub("_", ">", opt$mutation_type)

# read in the JSON that directs R/qtl2 to sample genotypes,
# phenotypes, and covariates
bxd <- read_cross2(opt$json)

Xcovar <- get_x_covar(bxd)

# insert pseudomarkers into the genotype map
# following lines are from https://kbroman.org/pages/teaching.html
gmap <- insert_pseudomarkers(bxd$gmap, step = 0.2, stepwidth = 'max')
pmap <- interp_map(gmap, bxd$gmap, bxd$pmap)

# read in the phenotype values for each BXD strain
df <- read.csv(opt$phenotype_file, header = TRUE)

# convert haplotype column to binary
phen_df <- df %>% mutate(haplotype_at_mutyh_binary = case_when(
    Haplotype_A == "D" ~ 1,
    Haplotype_A == "B" ~ 0,
))
phen_df <- phen_df %>% mutate(haplotype_at_ogg1_binary = case_when(
    Haplotype_B == "D" ~ 1,
    Haplotype_B == "B" ~ 0,
))

# calculate QTL genotype probabilities
pr <- calc_genoprob(bxd, pmap, error_prob = 0.002, map_function = "c-f")

# calculate kinship between strains using the
# "leave one chromosome out" method
k <- calc_kinship(pr, "loco")

# subset the phenotype data for this trial and convert to matrix
phen_df_sub <- subset(phen_df, mut == mutation_type & haplotype_at_mutyh_binary == 1)
phen_matrix <- as.matrix(phen_df_sub["Fraction"])

strain_names <- phen_df_sub$sample
rownames(phen_matrix) <- strain_names

covariate_matrix = as.matrix(phen_df_sub["haplotype_at_mutyh_binary"])
rownames(covariate_matrix) <- strain_names

# perform a genome scan, accounting for kinship and
# epoch as an additive covarirate
out <- scan1(pr, phen_matrix, kinship = k, Xcovar = Xcovar)

# perform a permutation test to assess significance
operm <- scan1perm(pr,
                   phen_matrix,
                   n_perm = 1000,
                   kinship = k,
                   Xcovar = Xcovar)

# get the LOD threshold for a < 0.05
alpha <- 0.05 / 7
lod_cutoff <- summary(operm, alpha = alpha)[1]

# print LOD peaks
peaks <- find_peaks(out, bxd$pmap, threshold = lod_cutoff)

# plot LOD scores genome-wide
png(opt$output_file, width = 7, height = 3.5, units = "in", res = 300)
par(mar = c(4.1, 4.1, 1.6, 1.1))
color <- c("black")
plot(out, pmap, lodcolumn = 1, col = color[1], bgcolor = "white", altbgcolor = "white", ylim = c(0, 6), main = mutation_type)
legend("topright", lwd = 1, col = "black", "LOD scores", bg = "white", lty = c(1, 1, 2))
abline(h = lod_cutoff, col = "darkgrey", lwd = 2, lty = 2)

dev.off()