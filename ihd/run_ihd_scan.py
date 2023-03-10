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
    compute_genotype_similarity,
    compute_manual_chisquare,
    compute_manual_cosine_distance,
)
from schema import IHDResultSchema, MutationSchema
import seaborn as sns
import statsmodels.api as sm
import numba


def main(args):

    # read in JSON file with file paths
    config_dict = None
    with open(args.config, "rb") as config:
        config_dict = json.load(config)

    # read in genotype info
    geno = pd.read_csv(config_dict['geno'])

    # read in singleton data and validate with pandera
    mutations = pd.read_csv(args.mutations)
    MutationSchema.validate(mutations)

    samples = mutations['sample'].unique()
    # get the overlap between those and the sample names in the genotype data
    samples_overlap = list(set(samples).intersection(set(geno.columns)))

    if len(samples_overlap) == 0:
        print("""Sorry, no samples in common between mutation data 
        and genotype matrix. Please ensure sample names are identical.""")
        sys.exit()

    # then subset the genotype and mutation data to include only those samples
    cols2use = ["marker"]
    cols2use.extend(samples_overlap)
    geno = geno[cols2use]

    mutations_filtered = mutations[mutations['sample'].isin(samples_overlap)]

    # for the null permutation test, shuffle the rows of the spectra
    # dataframe every time. otherwise keep it the same.
    samples, _, spectra = compute_spectra(mutations_filtered, k=args.k)
    print(f"""Using {len(samples)} samples and {int(np.sum(spectra))} 
    total mutations.""")

    # convert string genotypes to integers based on config definition
    replace_dict = config_dict['genotypes']
    geno_asint = geno.replace(replace_dict).replace({1: np.nan})

    # calculate allele frequencies at each site
    ac = np.nansum(geno_asint[samples], axis=1)
    an = np.sum(~np.isnan(geno_asint[samples]), axis=1) * 2
    afs = ac / an

    # only consider sites where allele frequency is between thresholds
    idxs2keep = np.where((afs > 0.1) & (afs < 0.9))[0]
    print("Using {} genotypes that meet filtering criteria.".format(idxs2keep.shape[0]))
    geno_filtered = geno_asint.iloc[idxs2keep][samples].values
    markers_filtered = geno_asint.iloc[idxs2keep]['marker'].values

    # compute similarity between allele frequencies at each marker 
    genotype_similarity = compute_genotype_similarity(geno_filtered)
    distance_method = compute_manual_chisquare
    if args.distance_method == "cosine":
        distance_method = compute_manual_cosine_distance

    # compute the maximum cosine distance between groups of
    # haplotypes at each site in the genotype matrix
    focal_dists = perform_ihd_scan(
        spectra,
        geno_filtered,
        genotype_similarity,
        distance_method=distance_method,
    )
    
    f, ax = plt.subplots()
    new_df = pd.DataFrame({
        "Genotype correlation": genotype_similarity,
        "Spectra distance": focal_dists,
    })
    sns.regplot(
        data=new_df,
        x="Genotype correlation",
        y="Spectra distance",
        ax=ax,
        scatter_kws={"alpha": 0.25}
    )
    X, y = new_df["Genotype correlation"].values, new_df["Spectra distance"].values
    lr = sm.OLS(y, sm.add_constant(X)).fit()
    #print (lr.summary())
    f.savefig("o.png", dpi=200)


    res_df = pd.DataFrame({
        'marker': markers_filtered,
        'Distance': focal_dists,
        'k': args.k,
    })
    IHDResultSchema.validate(res_df)

    # then do permutations
    null_distances = perform_permutation_test(
        spectra,
        geno_filtered,
        genotype_similarity,
        distance_method=distance_method,
        n_permutations=args.permutations,
        comparison_wide=args.comparison_wide,
    )

    # compute the 95th percentile of the maximum distance
    # distribution to figure out the distance corresponding to an alpha
    # of 0.05
    for pctile in (20, 5, 1):
        score_pctile = np.percentile(null_distances, 100 - pctile, axis=0)
        if not args.comparison_wide: score_pctile = score_pctile[0]
        res_df[f'{100 - pctile}th_percentile'] = score_pctile

    res_df.to_csv(args.out, index=False)

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
    p.add_argument(
        "-comparison_wide",
        action="store_true",
        help="Use per-marker instead of genome-wide distance thresholds.",
    )
    p.add_argument("-adj_column", default=None)
    p.add_argument("-chrom", default=None)
    p.add_argument("-distance_method", default="chisquare", type=str)
    args = p.parse_args()

    main(args)
