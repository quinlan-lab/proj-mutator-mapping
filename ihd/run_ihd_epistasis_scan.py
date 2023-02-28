import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import argparse
import json
import sys
from utils import (
    compute_spectra,
    perform_ihd_epistasis_scan,
    perform_epistasis_permutation_test,
)
from schema import IHDResultSchema, MutationSchema, MarkerMetadataSchema
from sklearn.manifold import MDS

def main(args):

    # read in JSON file with file paths
    config_dict = None
    with open(args.config, "rb") as config:
        config_dict = json.load(config)

    # read in genotype info
    geno = pd.read_csv(config_dict['geno'])

    # do the same with the genotype marker metadata
    markers = pd.read_csv(config_dict['markers'], dtype={"marker": str, "chromosome": str})
    MarkerMetadataSchema.validate(markers)
    geno = geno.merge(markers, on="marker")
    geno = geno[geno["chromosome"].isin(list(map(str, range(4, 8))))]

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

    # compute the maximum cosine distance between groups of
    # haplotypes at each site in the genotype matrix
    focal_dists = perform_ihd_epistasis_scan(
        spectra,
        geno_filtered,
    )



    dist_df = pd.DataFrame(focal_dists)
    dist_df["marker"] = markers_filtered
    dist_df = dist_df.merge(markers, on="marker")

    perm_dists = perform_epistasis_permutation_test(
        spectra,
        geno_filtered,
        n_permutations=args.permutations,
    )
    # get 95th pctile of dists
    pctile = np.percentile(perm_dists, 95)
    dist_df["threshold"] = pctile
    dist_df.to_csv(args.out, index=False)
    focal_dists[focal_dists <= pctile] = np.nan
    # pass_idxs = np.where(focal_dists >= pctile)[0][0]

    # pass_markers = markers_filtered.iloc[pass_idxs]

    # pass_markers["epi_dist"] =

    f, ax = plt.subplots()
    sns.heatmap(focal_dists, ax=ax)
    f.savefig("o.png", dpi=300)

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
