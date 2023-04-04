library(qtl2)
library(optparse)

option_list <- list(
  make_option(c("-j", "--json"), type="character", default=NULL),
  make_option(c("-p", "--phenotype_file"), type="character", default=NULL),
  make_option(c("-o", "--output_file"), type="character", default=NULL),
  make_option(c("-m", "--mutation_type", type="character", default=NULL)))

opt_parser <- OptionParser(option_list = option_list)
opt = parse_args(opt_parser)

mutation_type = gsub("_", ">", opt$mutation_type)

# read in the JSON that directs R/qtl2 to sample genotypes,
# phenotypes, and covariates
bxd <- read_cross2(opt$json)

# insert pseudomarkers into the genotype map
# following lines are from https://kbroman.org/pages/teaching.html
gmap <- insert_pseudomarkers(bxd$gmap, step=0.2, stepwidth='max')
pmap <- interp_map(gmap, bxd$gmap, bxd$pmap)

# read in the phenotype values for each BXD strain
phen_df <- read.csv(opt$phenotype_file, header = T)

# calculate QTL genotype probabilities
pr <- calc_genoprob(bxd, bxd$pmap, error_prob = 0.002, map_function = "c-f")

# calculate kinship between strains using the
# "leave one chromosome out" method
k = calc_kinship(pr, 'loco')

# subset the phenotype data for this trial and convert to matrix
phen_df_sub = subset(phen_df, Haplotype_A == "D" & mut == mutation_type)

phen_matrix <- as.matrix(phen_df_sub["CLR_fraction"])

strain_names <- phen_df_sub$sample
rownames(phen_matrix) <- strain_names

# perform a genome scan, accounting for kinship and
# epoch as an additive covarirate
out <- scan1(pr, phen_matrix)

# perform a permutation test to assess significance
operm <- scan1perm(pr, phen_matrix, n_perm=1000)

# get the LOD threshold for a < 0.05
alpha = 0.05 / 7
lod_cutoff = summary(operm, alpha=alpha)[1]

# print LOD peaks
peaks = find_peaks(out, bxd$pmap, threshold=lod_cutoff)

# plot LOD scores genome-wide
png(opt$output_file, width=7, height=3.5, units="in", res=300)
par(mar=c(4.1, 4.1, 1.6, 1.1))
color <- c("cornflowerblue")
plot(out, pmap, lodcolumn=1, col=color[1], ylim=c(0, 6), main=mutation_type)
legend("topright", lwd=1, col="cornflowerblue", "LOD scores", bg="gray90", lty=c(1,1,2))
abline(h=lod_cutoff, col='black', lwd=2, lty=2)

dev.off()