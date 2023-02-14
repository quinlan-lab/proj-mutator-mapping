import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import argparse
import json
import sys
from utils import (
    compute_spectra,
    perform_permutation_test,
    perform_ihd_scan,
)
from scipy.stats import ks_2samp, mannwhitneyu, ranksums
import seaborn as sns
from schema import IHDResultSchema, MutationSchema

plt.rc("font", size=18)

def main(args):

    # read in JSON file with file paths
    config_dict = None
    with open(args.config, "rb") as config:
        config_dict = json.load(config)

    # read in genotype info
    geno = pd.read_csv(config_dict['geno'])

    # read in marker info
    markers = pd.read_csv(config_dict["markers"])
    geno = geno.merge(markers, on="marker")

    # read in singleton data and validate with pandera
    mutations = pd.read_csv(args.mutations)
    MutationSchema.validate(mutations)

    samples = mutations['sample'].unique()
    # get the overlap between those and the sample names in the genotype data
    samples_overlap = list(set(samples).intersection(set(geno.columns)))

    if len(samples_overlap) == 0:
        print ("""Sorry, no samples in common between mutation data 
        and genotype matrix. Please ensure sample names are identical.""")
        sys.exit()

    # then subset the genotype and mutation data to include only those samples
    cols2use = ["marker", "chromosome"]
    cols2use.extend(samples_overlap)
    geno = geno[cols2use]

    mutations_filtered = mutations[mutations['sample'].isin(samples_overlap)]

    # for the null permutation test, shuffle the rows of the spectra
    # dataframe every time. otherwise keep it the same.
    samples, _, spectra = compute_spectra(mutations_filtered, k=args.k)

    print (f"""Using {len(samples)} samples and {int(np.sum(spectra))} 
    total mutations.""")

    # convert string genotypes to integers based on config definition
    replace_dict = config_dict['genotypes']

    res = []

    #geno['chromosome'] = "all"

    for chrom, chrom_geno in geno.groupby("chromosome"):

        geno_asint = chrom_geno.replace(replace_dict).replace({1: np.nan})

        # calculate allele frequencies at each site
        ac = np.nansum(geno_asint[samples], axis=1)
        an = np.sum(~np.isnan(geno_asint[samples]), axis=1) * 2
        afs = ac / an

        # only consider sites where allele frequency is between thresholds
        idxs2keep = np.where((afs > 0.1) & (afs < 0.9))[0]
        print ("Using {} genotypes that meet filtering criteria.".format(idxs2keep.shape[0]))
        geno_filtered = geno_asint.iloc[idxs2keep][samples].values
        markers_filtered = geno_asint.iloc[idxs2keep]['marker'].values
        chroms_filtered = geno_asint.iloc[idxs2keep]['chromosome'].values

        # compute the maximum cosine distance between groups of
        # haplotypes at each site in the genotype matrix
        focal_dists = perform_ihd_scan(
            spectra,
            geno_filtered,
        )

        # then do permutations
        max_distances = perform_permutation_test(
            spectra,
            geno_filtered,
            n_permutations=args.permutations,
            return_all=True,
        )

        # compute the k-s statistic between the focal and permuted distances
        ks_res = ks_2samp(focal_dists, max_distances, alternative="less")
        #ks_res = ranksums(focal_dists, max_distances, alternative="less", nan_policy="omit")
        res.append({
            "chromosome": chrom,
            "ks_stat": ks_res.statistic,
            "ks_p": -1 * np.log10(ks_res.pvalue),
        })

        f, ax = plt.subplots(figsize=(8, 6))
        max_val = max(focal_dists)
        if max(max_distances) > max_val: max_val = max(max_distances)
        bins=np.linspace(0, max_val, 20)
        colors = ["#E6803C", "#398D84"]
        labels = ["True", "Simulated"]
        ind = np.arange(bins[1:].shape[0])
        for i, d in enumerate((focal_dists, max_distances)):
            hist, edges = np.histogram(d, bins=bins)
            ax.bar(
                ind + (0.35 * i),
                hist / np.sum(hist),
                0.35,
                color=colors[i],
                label=labels[i],
                ec='k',
                lw=0.5,
            )
        ax.set_xticks(ind[::2])
        ax.set_xticklabels(ind[::2])
        ax.set_xlabel("Chi-square statistic bin (lower" + r"$\rightarrow$" + "higher)")
        ax.set_ylabel("Proportion of markers in bin")
        ax.set_title(f"Chi-square statistics between\nmarkers on chromosome {chrom}")
        ax.legend(frameon=False)
        f.tight_layout()
        f.savefig(f'csv/tmp/{chrom}.png', dpi=300)

    res_df = pd.DataFrame(res).query('chromosome != "X"').astype({"chromosome": int})
    f, (ax1, ax2) = plt.subplots(2, 1, figsize=(14, 16))
    sns.barplot(data=res_df, x="chromosome", y="ks_stat", ec='k', ax=ax1)
    sns.barplot(data=res_df, x="chromosome", y="ks_p", ec='k', ax=ax2)

    ax1.set_ylabel("K-S test statistic")
    ax2.set_ylabel("K-S test -log10(p)")

    ax2.set_xlabel("Chromosome")
    ax1.set_title("Comparisons of Chi-square statistic distributions between\ntrue and permuted mutation spectra on each chromosome")
    f.tight_layout()
    f.savefig("o.png", dpi=300)
    res_df.to_csv(args.out, index=False)

    f, ax = plt.subplots()


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument(
        "--mutations",
        type=str,
        help="Path to mutation data in CSV format.",
    )
    p.add_argument(
        "--config",
        type=str,
        help="Path to config file in JSON format.",
    )
    p.add_argument(
        "--out",
        help="Path in which to store the results of the IHD scan.",
    )
    p.add_argument(
        "-k",
        type=int,
        default=1,
        help="kmer size used for grouping mutations",
    )
    p.add_argument(
        "-permutations",
        type=int,
        default=1_000,
        help=
        "Number of permutations to perform when calculating significance thresholds.",
    )
    p.add_argument("-adj_column", default=None)
    p.add_argument("-chrom", default=None)
    args = p.parse_args()

    main(args)
