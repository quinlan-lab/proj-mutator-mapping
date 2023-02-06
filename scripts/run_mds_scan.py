import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import argparse
import json
import sys
from utils import (
    compute_spectra,
    compute_pairwise_marker_distance,
)
from sklearn.manifold import MDS
from skbio.stats.ordination import pcoa
from sklearn.metrics.pairwise import cosine_distances
from schema import IHDResultSchema, MutationSchema

def main(args):

    # read in JSON file with file paths
    config_dict = None
    with open(args.config, "rb") as config:
        config_dict = json.load(config)

    # read in genotype info
    geno = pd.read_csv(config_dict['geno'])
    markers = pd.read_csv(config_dict["markers"])

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
    cols2use = ["marker"]
    cols2use.extend(samples_overlap)
    geno = geno[cols2use]

    geno = geno.merge(markers, on="marker")

    mutations_filtered = mutations[mutations['sample'].isin(samples_overlap)]

    # for the null permutation test, shuffle the rows of the spectra
    # dataframe every time. otherwise keep it the same.
    samples, _, spectra = compute_spectra(mutations_filtered, k=args.k)

    print (f"""Using {len(samples)} samples and {int(np.sum(spectra))} 
    total mutations.""")
    res_df = []

    for chrom, geno_chrom in geno.groupby('chromosome'):

        # convert string genotypes to integers based on config definition
        replace_dict = config_dict['genotypes']
        geno_asint = geno_chrom.replace(replace_dict).replace({1: np.nan})

        # calculate allele frequencies at each site
        ac = np.nansum(geno_asint[samples], axis=1)
        an = np.sum(~np.isnan(geno_asint[samples]), axis=1) * 2
        afs = ac / an

        # only consider sites where allele frequency is between thresholds
        idxs2keep = np.where((afs > 0.1) & (afs < 0.9))[0]
        print ("Using {} genotypes that meet filtering criteria.".format(idxs2keep.shape[0]))
        geno_filtered = geno_asint.iloc[idxs2keep][samples].values
        markers_filtered = geno_asint.iloc[idxs2keep]['marker'].values

        # compute the pairwise_distances between the aggregate spectra
        # of haplotypes with the A or B (or any) allele at each pair of sites
        for genotype in (0, 2):
            spectra_sums = compute_pairwise_marker_distance(
                spectra,
                geno_filtered,
                genotype,
            )
            pairwise_dists = cosine_distances(spectra_sums)
            if chrom in ("4", "6"):
                f, ax = plt.subplots()
                sns.heatmap(pairwise_dists, ax=ax)
                f.savefig(f"{chrom}.{genotype}.heat.png", dpi=200)
            # then run MDS to get scaled distances at each marker
            n_comp = 10
            clf = MDS(n_components=n_comp, dissimilarity="precomputed", normalized_stress=False)
            mds_ = clf.fit_transform(pairwise_dists)
            #print (mds_)
            mds_df = pd.DataFrame(mds_, columns=[f"MDS{n}" for n in range(1, n_comp + 1)])
            mds_df["genotype"] = genotype
            mds_df['marker'] = markers_filtered
            res_df.append(mds_df)

    res_df = pd.concat(res_df)

    

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
    p.add_argument("-adj_column", default=None)
    p.add_argument("-chrom", default=None)
    args = p.parse_args()

    main(args)
