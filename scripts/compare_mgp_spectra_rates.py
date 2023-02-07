import argparse
import scipy.stats as ss
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from cyvcf2 import VCF
from statannotations.Annotator import Annotator

plt.rc('font', size=16)

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

smp2mutyh = dict(zip(mutyh_geno_formatted["strain"], mutyh_geno_formatted["Haplotype"]))

dumont = pd.read_excel(args.dumont_xls, sheet_name="TableS3", header=2)

mutations = [
         "C>A",
         "C>T",
         "C>G",
         "T>A",
         "T>C",
         "T>G"
     ]

cols2use = ["Strain"]
cols2use.extend([f"{m}.1" for m in mutations])

dumont_filtered = dumont[cols2use]

def make_rename_dict(colnames):
    revcomp = {"A": "T", "T": "A", "C": "G", "G": "C"}
    rename_dict = {}
    for c in colnames:
        orig, new = c.split('.')[0].split('>')
        if orig == "T":
            orig, new = "A", revcomp[new]
        new_c = r"$\rightarrow$".join([orig, new])
        rename_dict[c] = new_c
    return rename_dict


rename_dict = make_rename_dict(cols2use[1:])

# rename columns
dumont_filtered.rename(
    columns=rename_dict,
    inplace=True,
)

# format strain names to be more readable
dumont_filtered['strain'] = dumont_filtered['Strain'].apply(
    lambda x: x.replace('/', '_'))


dumont_merged = dumont_filtered.merge(sanger_genos, on="strain").drop(columns="Strain")

dumont_merged["Haplotypes at chr4 and chr6 peaks"] = dumont_merged["strain"].apply(lambda s: "-".join([smp2mutyh[s], smp2ogg[s]]))

dumont_merged = dumont_merged[dumont_merged["Haplotypes at chr4 and chr6 peaks"].isin([
    "B-B",
    "B-D",
    "D-B",
    "D-D",
])].drop(columns=["Haplotype"])


dumont_tidy = dumont_merged.melt(
    id_vars=["strain", "Haplotypes at chr4 and chr6 peaks"],
    var_name="Mutation type",
    value_name="Rate",
)


f, ax = plt.subplots(figsize=(10, 6))

palette = ["#398D84", "#E67F3A", "#EBBC2C", "#2F294A"]

res = []

sns.boxplot(
    data=dumont_tidy,
    x="Mutation type",
    y="Rate",
    hue="Haplotypes at chr4 and chr6 peaks",
    ax=ax,
    color="w",
    fliersize=0,
    linewidth=1,
)
sns.stripplot(
    data=dumont_tidy,
    x="Mutation type",
    y="Rate",
    hue="Haplotypes at chr4 and chr6 peaks",
    ax=ax,
    palette=palette,
    dodge=True,
    ec="k",
    linewidth=1,
)

ax.set_ylabel('Mutation rate (per bp, per gen.)')

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
    title="Haplotypes at chr4 and chr6 peaks",
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
    hue="Haplotypes at chr4 and chr6 peaks",
)
annotator.configure(
    test='t-test_welch',
    #test ="Mann-Whitney",
    comparisons_correction="BH",
    text_format='star',
)
#annotator.apply_and_annotate()

sns.set_style("ticks")
f.savefig(args.out + ".eps", dpi=300)
