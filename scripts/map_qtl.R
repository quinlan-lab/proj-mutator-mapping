library(qtl2)
library(optparse)

option_list <- list(
  make_option(c("-j", "--json"), type = "character", default = NULL),
  make_option(c("-p", "--phenotype_file"), type = "character", default = NULL),
  make_option(c("-o", "--output_file"), type = "character", default = NULL),
  make_option(c("-m", "--mutation_type"), type = "character", default = NULL),
  make_option(c("-t", "--trial"), type = "character", default = NULL))

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# read in the JSON that directs R/qtl2 to sample genotypes,
# phenotypes, and covariates
sim_bxd <- read_cross2(opt$json)

# read in the phenotype values for each BXD strain
phen_df <- read.csv(opt$phenotype_file, header = TRUE)

# calculate QTL genotype probabilities
pr <- calc_genoprob(sim_bxd, sim_bxd$pmap, error_prob = 0, map_function = "c-f")

# store results of simulated QTL scans
res_df <- data.frame(trial = -1,
                     bonferroni_corr = 1,
                     Power = 0)


# subset the phenotype data for this trial and convert to matrix
phenotype <- c(opt$mutation_type)
phen_matrix <- as.matrix(phen_df[phenotype])


strain_names <- phen_df$sample
rownames(phen_matrix) <- strain_names

# perform a genome scan, accounting for kinship and
# epoch as an additive covarirate
out <- scan1(pr, phen_matrix)

# perform a permutation test to assess significance
operm <- scan1perm(pr, phen_matrix, n_perm = 1000)

for (bc in c(1, 7, 96)) {

    # get the LOD threshold for a < 0.05
    alpha <- 0.05 / bc
    lod_cutoff <- summary(operm, alpha = alpha)[1]

    # print LOD peaks
    peaks <- find_peaks(out, sim_bxd$pmap, threshold = lod_cutoff)
    n_peaks <- dim(peaks)[1]
    discovered <- 0
    if (n_peaks > 0) {
        # markers at which a QTL was found
        predicted_markers <- as.vector(peaks[, 4])
        focal_marker_idx <- phen_df$focal_marker[1] + 1
        true_marker <- sim_bxd$pmap$`1`[focal_marker_idx]

        if (true_marker %in% predicted_markers) {
            discovered <- 1
        }

    }

    bc_df <- data.frame(trial = opt$trial,
                bonferroni_corr = bc,
                Power = discovered)

    res_df <- rbind(res_df, bc_df)
}



write.csv(subset(res_df, trial != -1), opt$output_file)
