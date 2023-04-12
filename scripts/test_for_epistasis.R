library(MASS)
library(ggplot2)
library(qtl2)
library(dplyr)
library(coxme)

df <- read.csv("/Users/tomsasani/quinlanlab/proj-mutator-mapping/csv/bxd_spectra_1mer.csv")

# get C>A mutations only
df <- subset(df, is_ca == 1)

# read in the JSON that directs R/qtl2 to sample genotypes,
# phenotypes, and covariates
bxd <- read_cross2("/Users/tomsasani/quinlanlab/proj-mutator-mapping/data/Rqtl_data/bxd.rqtl2.json")

# subset cross2 to relevant BXDs
bxd <- bxd[df$sample, ]

# calculate QTL genotype probabilities
# (error_prob taken directly from R/qtl2 user guide)
pr <- calc_genoprob(bxd, bxd$pmap, error_prob = 0.002, map_function = "c-f")

# calculate kinship between strains using the
# "overall" method
k <- calc_kinship(pr, "overall")

# convert haplotype column to binary
phen_df <- df %>% mutate(haplotype_at_mutyh_binary = case_when(
    Haplotype_A == "D" ~ 1,
    Haplotype_A == "B" ~ 0,
))
phen_df <- phen_df %>% mutate(haplotype_at_ogg1_binary = case_when(
    Haplotype_B == "D" ~ 1,
    Haplotype_B == "B" ~ 0,
))

kinship <- as.matrix(k)
colnames(kinship) <- phen_df$sample
rownames(kinship) <- phen_df$sample

# https://sahirbhatnagar.com/blog/2017/10/12/mixed-models-with-kinship-in-r/

# run lmekin to test for interaction effects between genotypes on
# C>A mutation fractions
gfit1 <- lmekin(CLR_fraction ~ haplotype_at_mutyh_binary + haplotype_at_ogg1_binary + (1 | sample),
    data = phen_df,
    varlist = kinship
)
gfit2 <- lmekin(CLR_fraction ~ haplotype_at_mutyh_binary * haplotype_at_ogg1_binary + (1 | sample),
    data = phen_df,
    varlist = kinship
)
print (gfit2)

# use a poisson model to test for interaction effects between genotypes
# on C>A mutation counts, modeled as rates
m1 <- glm(Count ~ offset(log(ADJ_AGE)) + haplotype_at_mutyh_binary + haplotype_at_ogg1_binary,
    data = phen_df,
    family = poisson()
)
m2 <- glm(Count ~ offset(log(ADJ_AGE)) + haplotype_at_mutyh_binary * haplotype_at_ogg1_binary,
    data = phen_df,
    family = poisson()
)
print(anova(m2, m1, test = "Chisq"))

# use the MGP data from dumont to test for interaction effects between
# genotypes on C>A mutation counts, modeled as rates
df <- read.csv("/Users/tomsasani/quinlanlab/proj-mutator-mapping/csv/mgp_spectra.csv")

df <- subset(df, is_ca == 1)
m1 <- glm(Count ~ offset(log(CALLABLE_C)) + Haplotype_A + Haplotype_B, 
    data = df, 
    family = poisson()
)
m2 <- glm(Count ~ offset(log(CALLABLE_C)) + Haplotype_A * Haplotype_B, 
    data = df, 
    family = poisson()
)

print(anova(m2, m1, test = "Chisq"))
