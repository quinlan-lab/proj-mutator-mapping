library(qtl2)

# read in the JSON that directs R/qtl2 to sample genotypes,
# phenotypes, and covariates
bxd <- read_cross2("/Users/tomsasani/quinlanlab/proj-mutator-mapping/data/Rqtl_data/bxd.rqtl2.json")

# read in the phenotype values for each BXD strain
phen_df = read.csv("/Users/tomsasani/quinlanlab/proj-mutator-mapping/csv/tidy_spectra.csv", header=T)
phen_df = subset(phen_df, sample != "BXD68")

# get phenotypes as a matrix 
phenotypes = c("C_A", "C_T", "C_G", "A_T", "A_C", "A_G")
phen_matrix = as.matrix(phen_df[phenotypes])


strain_names = phen_df$sample
rownames(phen_matrix) = strain_names

# subset cross2 to relevant BXDs
bxd = bxd[phen_df$sample, ]

# insert pseudomarkers into the genotype map
# following lines are from https://kbroman.org/pages/teaching.html
gmap <- insert_pseudomarkers(bxd$gmap, step=0.2, stepwidth='max')
pmap <- interp_map(gmap, bxd$gmap, bxd$pmap)

# calculate QTL genotype probabilities
# (error_prob taken directly from R/qtl2 user guide)
pr <- calc_genoprob(bxd, gmap, error_prob=0.002, map_function="c-f")

# get covariates as a matrix 
covariate_cols = c("haplotype_at_chr4_qtl")#, "haplotype_at_chr6_qtl")
covariate_matrix = as.matrix(phen_df[covariate_cols])
rownames(covariate_matrix) = phen_df$sample

res_df <- NULL

# loop over chromosomes and calculate heritability 
for (chrom in 1:19) {

    chrom_str = toString(chrom)

    # calculate kinship between strains using the
    # "overall" method
    k = calc_kinship(pr[, chrom], 'overall')

    # perform genome scan
    hsq <- est_herit(phen_matrix, k)#, covariate_matrix)

    print (chrom)
    print (hsq)

    idxs = c(1:6)
        
    if (is.null(res_df)) {
        res_df <- data.frame(chromosome=rep(chrom_str, 6),
                        heritability=hsq[idxs],
                        mutation=phenotypes)
        
    }
    else {
        hsq_df <- data.frame(chromosome=rep(chrom_str, 6),
                        heritability=hsq[idxs],
                        mutation=phenotypes)

        res_df <- rbind(res_df, hsq_df)
    }

    
}

write.csv(res_df, "/Users/tomsasani/quinlanlab/proj-mutator-mapping/csv/heritability.csv")
