import argparse
import scipy.stats as ss
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from cyvcf2 import VCF
from statannotations.Annotator import Annotator

plt.rc('font', size=20)

p = argparse.ArgumentParser()
p.add_argument(
    "--dumont_xls",
    required=True,
    help="""strain-private substitution data from Dumont et al. (2019) \
                        as it appears in Table 3 of the Supplementary Data""",
)
p.add_argument("--vcf")
p.add_argument(
    "--out",
    required=True,
    help="""name of output plot""",
)
p.add_argument("--mutyh")

args = p.parse_args()

# figure out which strains possess the Ogg1 mutation
mutation_site = f"6:113328509-113328510"
mutation_site = f"6:115849643-115849644"
vcf = VCF(args.vcf, gts012=True)
sanger_genos = []

for v in vcf(mutation_site):
    dp4 = v.format("DP4")
    ad, rd = np.sum(dp4[:, :2], axis=1), np.sum(dp4[:, 2:], axis=1)
    td = ad + rd
    gts = v.gt_types
    for si, s in enumerate(vcf.samples):
        gt = -1
        if td[si] >= 10:
            gt = gts[si]
        sanger_genos.append({
            "strain": s,
            "Haplotype": "B" if gt == 0 else "D" if gt == 2 else "UNK"
        })

sanger_genos = pd.DataFrame(sanger_genos).query('Haplotype != "UNK"')
smp2ogg = dict(zip(sanger_genos["strain"], sanger_genos["Haplotype"]))
mutyh_genos = pd.read_csv(args.mutyh, names=["site", "strain", "genotype", "ref", "alt"], dtype={"genotype": str})

genotype2hap = {
    "0,0,0,0,0": "B",
    "2,0,0,2,2": "intermediate",
    "2,2,2,2,2": "D",
}

mutyh_geno_formatted = mutyh_genos.groupby("strain").agg({"genotype": lambda s: genotype2hap[",".join(s)] if ",".join(s) in genotype2hap else "UNK"}).rename(columns={"genotype": "Haplotype"}).reset_index()

smp2mutyh = dict(
    zip(
        mutyh_geno_formatted["strain"],
        mutyh_geno_formatted["Haplotype"],
    ))

# read in numbers of callable base pairs in each strain
dumont_callable = pd.read_excel(
    args.dumont_xls,
    sheet_name="TableS3",
    header=2,
)

cols2use = ["Strain"]
cols2use.extend(["nC", "nT", "nA", "nG"])

dumont_callable = dumont_callable[cols2use]
dumont_callable['strain'] = dumont_callable['Strain'].apply(
    lambda x: x.replace('/', '_'))


# read in per-sample mutation data
dumont_long = pd.read_excel(args.dumont_xls, sheet_name="TableS4", header=2)
revcomp = {"A": "T", "T": "A", "C": "G", "G": "C"}

def convert_to_mutation(row):
    orig, new = row["ancestral"], row["derived"]
    fp, tp = row["5prime_flank"], row["3prime_flank"]
    mutation = None
    # if we need to reverse complement
    if orig in ("G", "T"):
        orig, new = revcomp[orig], revcomp[new]
        fp, tp = revcomp[fp], revcomp[tp]
    if orig == "C" and tp == "G":
        mutation = r"$\rightarrow$".join(["CpG", "TpG"])
        #mutation = r"$\rightarrow$".join(["C", "T"])
    else:
        mutation = r"$\rightarrow$".join([orig, new])
    if mutation is None: print (row["FocalStrain"], row["chr"], row["pos"])
    return mutation

# add formatted mutation column to each row
dumont_long["Mutation type"] = dumont_long.apply(lambda row: convert_to_mutation(row), axis=1)

dumont_long_grouped = dumont_long.groupby(["FocalStrain", "Mutation type"]).size().reset_index().rename(columns={"FocalStrain": "strain", 0: "Count"})
# merge per-sample mutation counts with callable base pairs
dumont_merged = dumont_long_grouped.merge(dumont_callable, on="strain")

# merge mutation data with genotypes at markers
dumont_merged = dumont_merged.merge(sanger_genos, on="strain").drop(columns="Strain")
dumont_merged["Haplotypes"] = dumont_merged["strain"].apply(lambda s: "-".join([smp2mutyh[s], smp2ogg[s]]))

dumont_merged["CALLABLE_C"] = dumont_merged["nC"] + dumont_merged["nG"]
dumont_merged["CALLABLE_A"] = dumont_merged["nA"] + dumont_merged["nT"]


dumont_merged.drop(columns=["nC", "nT", "nA", "nG", "Haplotype"], inplace=True)


dumont_merged["Rate"] = dumont_merged.apply(
    lambda row: row["Count"] / row["CALLABLE_C"] if row["Mutation type"].
    startswith("C") else row["Count"] / row["CALLABLE_A"],
    axis=1,
)

# sum the adjusted per-bp rates in each strain
total_rates = dumont_merged.groupby(["strain"]).agg({"Rate": sum}).reset_index().rename(columns={"Rate": "Total rate"})
dumont_tidy = dumont_merged.merge(total_rates, on="strain")

dumont_tidy["Fraction"] = dumont_tidy["Rate"] / dumont_tidy["Total rate"]
print (dumont_tidy)
dumont_tidy = dumont_tidy[dumont_tidy["Haplotypes"].isin(["B-B", "B-D", "D-B", "D-D"])]
dumont_tidy["is_ca"] = dumont_tidy["Mutation type"].apply(lambda m: 1 if m == r"C$\rightarrow$A" else 0)
dumont_tidy["Haplotype_A"] = dumont_tidy["Haplotypes"].apply(lambda h: h.split('-')[0])
dumont_tidy["Haplotype_B"] = dumont_tidy["Haplotypes"].apply(lambda h: h.split('-')[1])
dumont_tidy.to_csv("/Users/tomsasani/quinlanlab/proj-mutator-mapping/csv/dumont_tidy.csv")

f, ax = plt.subplots(figsize=(14, 6))

palette = dict(
        zip(
            ["B-B", "B-D", "D-B", "D-D"],
            ["#398D84", "#E67F3A", "#EBBC2C", "#2F294A"],
        ))
res = []

sns.boxplot(
    data=dumont_tidy.sort_values(["Mutation type", "Haplotypes"], ascending=True),
    x="Mutation type",
    y="Fraction",
    hue="Haplotypes",
    ax=ax,
    color="w",
    fliersize=0,
    linewidth=1,
)
sns.stripplot(
    data=dumont_tidy.sort_values(["Mutation type", "Haplotypes"], ascending=True),
    x="Mutation type",
    y="Fraction",
    hue="Haplotypes",
    ax=ax,
    palette=palette,
    dodge=True,
    ec="k",
    linewidth=1,
)

ax.set_ylabel('')#'Mutation rate (per bp, per gen.)')
ax.set_ylim(0, 0.6)
# change all spines
# for axis in ['top','bottom','left','right']:
#     ax.spines[axis].set_linewidth(2.5)

# increase tick width
#ax.tick_params(width=2.5)
sns.despine(ax=ax, top=True, right=True)

handles, labels = ax.get_legend_handles_labels()
l = plt.legend(
    handles[4:],
    labels[4:],
    title="Genotypes at chr4 and chr6 peaks",
    frameon=False,
)


res = []

mutation_type = r"$\rightarrow$".join(["C", "A"])
pairs = [
    ((mutation_type, "B-B"), (mutation_type, "B-D")),
    ((mutation_type, "B-D"), (mutation_type, "D-D")),
    ((mutation_type, "D-B"), (mutation_type, "D-D")),
    ((mutation_type, "B-B"), (mutation_type, "D-B")),
    ((mutation_type, "B-B"), (mutation_type, "D-D")),

]

annotator = Annotator(
    ax,
    pairs,
    data=dumont_tidy,
    x="Mutation type",
    y="Rate",
    hue="Haplotypes",
)
annotator.configure(
    test='t-test_welch',
    #test ="Mann-Whitney",
    comparisons_correction="BH",
    text_format='star',
)
#annotator.apply_and_annotate()
# change all spines
for axis in ['top','bottom','left','right']:
    ax.spines[axis].set_linewidth(1.5)

# increase tick width
ax.tick_params(width=1.5)
#sns.set_style("ticks")
f.tight_layout()
f.savefig(args.out + ".png", dpi=300)
f.savefig(args.out + ".eps")