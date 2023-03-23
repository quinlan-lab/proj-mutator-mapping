import argparse
import scipy.stats as ss
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from typing import Dict
import seaborn as sns
from cyvcf2 import VCF
from statannotations.Annotator import Annotator

REV_COMP = {
    "A": "T",
    "T": "A",
    "C": "G",
    "G": "C",
}


def convert_to_mutation(row):
    orig, new = row["ancestral"], row["derived"]
    fp, tp = row["5prime_flank"], row["3prime_flank"]
    mutation = None
    # if we need to reverse complement
    if orig in ("G", "T"):
        orig, new = REV_COMP[orig], REV_COMP[new]
        fp, tp = REV_COMP[fp], REV_COMP[tp]
    if orig == "C" and tp == "G":
        mutation = r"$\rightarrow$".join(["CpG", "TpG"])
        #mutation = r"$\rightarrow$".join(["C", "T"])
    else:
        mutation = r"$\rightarrow$".join([orig, new])
    if mutation is None: print (row["FocalStrain"], row["chr"], row["pos"])
    return mutation

def get_genotype_at_site(vcf: VCF, pos: str) -> Dict[str, str]:
    genos = []
    for v in vcf(pos):
        dp4 = v.format("DP4")
        ad, rd = np.sum(dp4[:, :2], axis=1), np.sum(dp4[:, 2:], axis=1)
        td = ad + rd
        gts = v.gt_types
        for si, s in enumerate(vcf.samples):
            gt = -1
            if td[si] >= 10:
                gt = gts[si]
            genos.append({
                "strain": s,
                "Haplotype": "B" if gt == 0 else "D" if gt == 2 else "UNK"
            })

    # remove samples with UNK genotypes at the site and create a dictionary
    # mapping samples to their genotypes at the locus
    geno_df = pd.DataFrame(genos).query('Haplotype != "UNK"')
    smp2geno = dict(zip(geno_df["strain"], geno_df["Haplotype"]))

    return smp2geno

plt.rc('font', size=20)

def main(args):

    # read in MGP VCF
    vcf = VCF(args.vcf, gts012=True)

    # figure out which MGP strains possess D- or B-like alleles
    # at the markers on chromosome 4 and chromosome 6
    chr4_marker_pos = "4:116750874-116750875"
    chr6_marker_pos = "6:114052413-114052414"

    chr4_genos = get_genotype_at_site(vcf, chr4_marker_pos)
    chr6_genos = get_genotype_at_site(vcf, chr6_marker_pos)

    # read in numbers of callable base pairs in each strain
    dumont_callable = pd.read_excel(
        args.dumont_xls,
        sheet_name="TableS3",
        header=2,
    )

    cols2use = ["Strain"]
    cols2use.extend(["nC", "nT", "nA", "nG"])

    # reformat column names
    dumont_callable = dumont_callable[cols2use]
    dumont_callable['strain'] = dumont_callable['Strain'].apply(
        lambda x: x.replace('/', '_'))

    # read in per-sample mutation data in MGP strains
    dumont_long = pd.read_excel(args.dumont_xls, sheet_name="TableS4", header=2)
    # add formatted mutation column to each row
    dumont_long["Mutation type"] = dumont_long.apply(lambda row: convert_to_mutation(row), axis=1)
    # figure out total counts of each mutation type in each MGP strain
    dumont_long_grouped = dumont_long.groupby([
        "FocalStrain",
        "Mutation type",
    ]).size().reset_index().rename(columns={
        "FocalStrain": "strain",
        0: "Count",
    })

    # merge per-sample mutation counts with callable base pairs
    dumont_merged = dumont_long_grouped.merge(dumont_callable, on="strain")
    # annotate merged dataframe with information about haplotypes at each marker
    dumont_merged["Haplotypes"] = dumont_merged["strain"].apply(lambda s: "-".join([chr4_genos[s], chr6_genos[s]]))
    # figure out total numbers of callable cytosines and adenines in each strain
    dumont_merged["CALLABLE_C"] = dumont_merged["nC"] + dumont_merged["nG"]
    dumont_merged["CALLABLE_A"] = dumont_merged["nA"] + dumont_merged["nT"]
    # and drop superfluous columns
    dumont_merged.drop(columns=["nC", "nT", "nA", "nG"], inplace=True)
    print (dumont_merged.groupby("Haplotypes").size())

    # convert mutation counts into rates, taking into account the nucleotide context
    # in which each 1-mer occurred
    dumont_merged["Rate"] = dumont_merged.apply(
        lambda row: row["Count"] / row["CALLABLE_C"] if row["Mutation type"].
        startswith("C") else row["Count"] / row["CALLABLE_A"],
        axis=1,
    )
    # sum the context-adjusted rates in each strain
    total_rates = dumont_merged.groupby(["strain"]).agg({"Rate": sum}).reset_index().rename(columns={"Rate": "Total rate"})

    # create a tidy dataframe with 1-mer and total rates in each strain
    dumont_tidy = dumont_merged.merge(total_rates, on="strain")
    # convert rates to fractions as well
    dumont_tidy["Fraction"] = dumont_tidy["Rate"] / dumont_tidy["Total rate"]
    # filter the tidy dataframe to only include strains with non-UNK genotypes
    # at both of the mutator loci
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

    # change all spines
    for axis in ['top','bottom','left','right']:
        ax.spines[axis].set_linewidth(1.5)

    # increase tick width
    ax.tick_params(width=1.5)
    f.tight_layout()
    f.savefig(args.out, dpi=300)

if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument(
        "--dumont_xls",
        help="""strain-private substitution data from Dumont et al. (2019) \
                            as it appears in Table 3 of the Supplementary Data""",
    )
    p.add_argument(
        "--vcf",
        help="VCF file containing relevant mutation data in MGP strains",
    )
    p.add_argument(
        "--out",
        help="""name of output plot""",
    )
    args = p.parse_args()
    main(args)