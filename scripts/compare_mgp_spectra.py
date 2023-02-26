import argparse
import scipy.stats as ss
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from cyvcf2 import VCF
from statannotations.Annotator import Annotator
from statsmodels.stats.multitest import multipletests

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
cols2use.extend([f"{m}" for m in mutations])
cols2use.extend(["nA", "nC", "nT", "nG"])

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


rename_dict = make_rename_dict(cols2use[1:7])

# rename columns
dumont_filtered.rename(
    columns=rename_dict,
    inplace=True,
)

# format strain names to be more readable
dumont_filtered['strain'] = dumont_filtered['Strain'].apply(
    lambda x: x.replace('/', '_'))

# map samples to corresponding indices
smp2idx = dict(
    zip(
        dumont_filtered['strain'],
        range(len(dumont_filtered['strain'])),
    ))

# get the counts of singletons of each mutation type
# in each strain
denoms = dumont_filtered.values[:, 7:]

# get the counts of callable base pairs of each nucleotide
# type in each strain
counts = dumont_filtered.values[:, 1:7]


dumont_merged = dumont_filtered.merge(sanger_genos, on="strain").drop(columns="Strain")

dumont_merged["Combined haplotype"] = dumont_merged["strain"].apply(lambda s: "-".join([smp2mutyh[s], smp2ogg[s]]))

dumont_merged = dumont_merged[dumont_merged["Combined haplotype"].isin([
    "B-B",
    "B-D",
    "D-B",
    "D-D",
])]

mutations = dumont_filtered.columns[1:7]
mut2idx = dict(zip(mutations, range(len(mutations))))
strain2sums = dict(zip(dumont_merged['strain'], np.sum(dumont_merged[mutations], axis=1)))

dumont_merged.drop(columns=["nA", "nT", "nC", "nG", "Haplotype"], inplace=True)

dumont_tidy = dumont_merged.melt(id_vars=["strain", "Combined haplotype"], var_name="Mutation type", value_name="Count")
dumont_tidy["Total"] = dumont_tidy["strain"].apply(lambda s: strain2sums[s])

dumont_grouped = dumont_tidy.groupby(["Combined haplotype", "Mutation type"]).agg({"Count": sum, "Total": sum}).reset_index()

dumont_grouped["Aggregate fraction"] = dumont_grouped["Count"] / dumont_grouped["Total"]


pairs, pvalues = [], []

# compare spectra between strains with different
# configurations of Mutyh mutations
for cat_a, cat_b in [
    ("D-D", "B-B"),
    #("D-D", "B-D"),
    #("D-B", "B-B"),
    ("B-D", "B-B"),
    ("D-D", "D-B")
]:

    a_smps = dumont_tidy[dumont_tidy["Combined haplotype"] == cat_a]["strain"].unique()
    b_smps = dumont_tidy[dumont_tidy["Combined haplotype"] == cat_b]["strain"].unique()
    # sample indices in each category
    a_idx = np.array([smp2idx[s] for s in a_smps if s in smp2idx])
    b_idx = np.array([smp2idx[s] for s in b_smps if s in smp2idx])

    # subsets of mutation counts for strains in each
    # mutation "category"
    subset_0 = counts[a_idx, :]
    subset_1 = counts[b_idx, :]

    # get the callable number of basepairs for each
    # nucleotide for the samples in either category
    subset_0_denom = denoms[a_idx, :]
    subset_1_denom = denoms[b_idx, :]

    # sum the mutation counts of all types for strains
    # in each category
    subset_0_totals = np.sum(subset_0, axis=0)
    subset_1_totals = np.sum(subset_1, axis=0)

    # get the total number of callable base pairs (across all
    # nucleotides) in either subset
    subset_0_denom_totals = np.sum(subset_0_denom, axis=0)
    subset_1_denom_totals = np.sum(subset_1_denom, axis=0)

    nuc2denom = {'A': 0, 'C': 1, 'G': 2, 'T': 3}

    nuc2comp = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}

    a_adj, b_adj = [], []

    # we want to perform a Chi-square test to compare mutation
    # spectra in subset 0 vs. subset 1, but we need to first
    # adjust the counts of each mutation type in each subset by
    # the number of A, C, T, or G base pairs that met filtering
    # criteria across all strains in either subset.
    for mut, mut_i in mut2idx.items():
        base_nuc = mut[0]

        # first, sum the total number of mutations of type `mut`
        # observed in each subset
        a_fore = subset_0_totals[mut_i]
        b_fore = subset_1_totals[mut_i]

        # then get thethe sum of callable base pairs that are
        # X nucleotides across strains in subset 0, and same for subset 1.
        # we add the sum of callable base pairs that
        # are the complement of X, as well
        a_denom = subset_0_denom_totals[nuc2denom[base_nuc]]
        a_denom += subset_0_denom_totals[nuc2denom[nuc2comp[base_nuc]]]
        b_denom = subset_1_denom_totals[nuc2denom[base_nuc]]
        b_denom += subset_1_denom_totals[nuc2denom[nuc2comp[base_nuc]]]

        # we now adjust the counts of each mutation type in each strain
        if a_denom > b_denom:
            adj = b_denom / a_denom
            a_fore = a_fore * adj
        elif b_denom > a_denom:
            adj = a_denom / b_denom
            b_fore = b_fore * adj

        a_adj.append(a_fore)
        b_adj.append(b_fore)

    # for every mutation type, perform a chi-square test of
    # enrichment of that mutation type's counts in subset 0 vs. subset 1
    for mut, mut_i in mut2idx.items():
        a_fore = a_adj[mut_i]
        a_back = sum(a_adj) - a_fore
        b_fore = b_adj[mut_i]
        b_back = sum(b_adj) - b_fore

        _, p, _, _ = ss.chi2_contingency([[a_fore, b_fore], [a_back, b_back]])
        if p >= 0.05: continue
        pairs.append(((mut, cat_a), (mut, cat_b)))
        pvalues.append(p)
        
palette = ["#398D84", "#E67F3A", "#EBBC2C", "#2F294A"]

res = []

plt.figure(figsize=(8, 6))

ax = sns.barplot(
    data=dumont_grouped,
    x="Mutation type",
    y="Aggregate fraction",
    hue="Combined haplotype",
    #ax=ax,
    palette=palette,
    ec="k",
    linewidth=2,
)

annotator = Annotator(
        ax,
        pairs,
        data=dumont_grouped,
        x="Mutation type",
        y="Aggregate fraction",
        hue="Combined haplotype",
    )


annotator.configure(
    test=None,
    test_short_name=r"$\chi^{2}$",
    #text_format="full",
    #text_format="simple",
    text_format="star",
    #comparisons_correction="BH",
).set_pvalues(pvalues=pvalues).annotate()

# change all spines
for axis in ['top','bottom','left','right']:
    ax.spines[axis].set_linewidth(1.5)

# increase tick width
ax.tick_params(width=1.5)

handles, labels = ax.get_legend_handles_labels()
l = plt.legend(
    handles,
    labels,
    title="Haplotypes at chr4 and chr6 peaks",
    frameon=False,
)

#sns.set_style("ticks")
sns.despine(ax=ax, top=True, right=True)

plt.tight_layout()
plt.savefig(args.out + ".eps")
plt.savefig(args.out + ".png", dpi=300)

