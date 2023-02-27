library(MASS)
library(ggplot2)

df <- read.csv("/Users/tomsasani/quinlanlab/proj-mutator-mapping/csv/true_tidy_spectra.csv")

df <- subset(df, is_ca == 1)


g <- ggplot(df, aes(x = Haplotype, y = BC_Fraction)) +
    geom_boxplot(outlier.size=0) +
    geom_jitter(aes(col = Haplotype), width = 0.2) +
    facet_wrap(~Epoch)

ggsave("o.png", g)

# df = read.csv("/Users/tomsasani/quinlanlab/proj-mutator-mapping/csv/dumont_tidy.csv")

# df = subset(df, is_ca == 1)

# m2 = lm(Rate ~ Haplotype_A + Haplotype_B, data=df)
# m3 = lm(Rate ~ Haplotype_A * Haplotype_B, data=df)
# print (anova(m2, m3))#, test="Chisq"))
# print (summary(m3))

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


cols2use <- c("haplotype_at_mutyh_binary", "haplotype_at_ogg1_binary", "Epoch", "sample")

x1 <- cbind(1, as.matrix(phen_df[cols2use]))

Y <- phen_df$CLR_Fraction

dat <- data.frame(Y, mutyh_gt = x1[, 2], ogg1_gt = x1[, 3], epoch = x1[, 4], sample = x1[, 5])

kinship <- as.matrix(k)
diag(k) <- 1

colnames(kinship) <- phen_df$sample
rownames(kinship) <- phen_df$sample

# https://sahirbhatnagar.com/blog/2017/10/12/mixed-models-with-kinship-in-r/

gfit1 <- lm(Y ~ mutyh_gt * ogg1_gt, data = dat)
gfit2 <- lmekin(Y ~ mutyh_gt * ogg1_gt + (1 | sample), data = dat, varlist = kinship)

print(summary(gfit1))
print(gfit2)

m1 <- glm(Count ~ offset(log(ADJ_AGE)) + Haplotype_A + Haplotype_B, data = df, family=poisson(link="log"))
m2 <- glm(Count ~ offset(log(ADJ_AGE)) + Haplotype_A * Haplotype_B, data = df, family=poisson(link="log"))
print (anova(m2, m1, test="Chisq"))
