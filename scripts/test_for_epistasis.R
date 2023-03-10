library(MASS)
library(ggplot2)

df <- read.csv("/Users/tomsasani/quinlanlab/proj-mutator-mapping/csv/true_tidy_spectra.csv")

df <- subset(df, is_ca == 1)

library(qtl2)
library(dplyr)
library(coxme)

# read in the JSON that directs R/qtl2 to sample genotypes,
# phenotypes, and covariates
bxd <- read_cross2("/Users/tomsasani/quinlanlab/proj-mutator-mapping/data/Rqtl_data/bxd.rqtl2.json")

# subset cross2 to relevant BXDs
bxd <- bxd[df$sample, ]

# insert pseudomarkers into the genotype map
# following lines are from https://kbroman.org/pages/teaching.html
gmap <- insert_pseudomarkers(bxd$gmap, step = 0.2, stepwidth = "max")

# calculate QTL genotype probabilities
# (error_prob taken directly from R/qtl2 user guide)
pr <- calc_genoprob(bxd, gmap, error_prob = 0.002, map_function = "c-f")

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

gfit1 <- lm(Fraction ~ Epoch + Haplotype_A * Haplotype_B, data = phen_df)
gfit1 <- lmekin(Fraction ~ Epoch + Haplotype_A + Haplotype_B + (1 | sample),
                 data = phen_df,
                 varlist = kinship)
gfit2 <- lmekin(Fraction ~ Haplotype_A * Haplotype_B + (1 | sample),
                 data = phen_df,
                 varlist = kinship)

print(gfit2)
#print (anova(gfit2, gfit1))

m1 <- glm(Count ~ offset(log(ADJ_AGE)) + Haplotype_A + Haplotype_B,
            data = df,
            family = poisson(link = "log"))
m2 <- glm(Count ~ offset(log(ADJ_AGE)) + Haplotype_A * Haplotype_B,
            data = df,
            family = poisson(link = "log"))
print (anova(m2, m1, test = "Chisq"))

df = read.csv("/Users/tomsasani/quinlanlab/proj-mutator-mapping/csv/dumont_tidy.csv")

df = subset(df, is_ca == 1)
m1 <- glm(Count ~ offset(log(CALLABLE_C)) + Haplotype_A + Haplotype_B, data = df, family=poisson(link="log"))
m2 <- glm(Count ~ offset(log(CALLABLE_C)) + Haplotype_A * Haplotype_B, data = df, family=poisson(link="log"))
print (anova(m2, m1, test="Chisq"))
