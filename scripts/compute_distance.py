import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import itertools
import argparse
import json
import sys
from utils import (
    compute_spectra,
    compute_haplotype_distance,
    perform_permutation_test,
    perform_ihd_scan,
)

def main(args):

    # read in JSON file with file paths
    config_dict = None
    with open(args.config, "rb") as config:
        config_dict = json.load(config)

    geno_raw = pd.read_csv(config_dict['geno'])
    marker_info = pd.read_csv(config_dict['markers'])

    # merge genotype and marker information
    geno = geno_raw.merge(marker_info, on="marker")
    geno = geno[geno['chromosome'] != "X"]


    singletons = pd.read_csv(args.singletons)
    singletons['count'] = 1

    if args.chrom:
        base_chrom = args.chrom
        if 'chr' in args.chrom:
            base_chrom = args.chrom[3:]

        geno = geno[(geno['chromosome'] == base_chrom) | (geno['chromosome'] == f"chr{base_chrom}") ]
        singletons = singletons[(singletons['chrom'] == base_chrom) | (singletons['chrom'] == f"chr{base_chrom}")]

    samples = singletons['Strain'].unique()
    # get the overlap between those and the sample names in the genotype data
    samples_overlap = list(set(samples).intersection(set(geno.columns)))

    if len(samples_overlap) == 0:
        print ("No samples in common between mutation data and genotype matrix.")
        sys.exit()

    # then subset the genotype data to include only those samples
    cols2use = ["marker"]
    cols2use.extend(samples_overlap)
    geno = geno[cols2use]

    singletons = singletons[singletons['Strain'].isin(samples_overlap)]

    # for the null permutation test, shuffle the rows of the spectra
    # dataframe every time. otherwise keep it the same.
    samples, mutations, spectra = compute_spectra(singletons, k=args.k)
    smp2idx = dict(zip(samples, range(len(samples))))

    # get sums of mutation counts in each strain
    spectra_sums = np.sum(spectra, axis=1)
    spectra2keep = np.where(spectra_sums >= 0)[0]
    samples, spectra = list(np.array(samples)[spectra2keep]), spectra[spectra2keep]

    print (f"Using {len(samples)} samples and {int(np.sum(spectra))} total mutations.")

    # weight spectra if desired
    weight_dict = None
    if args.adj_column:
        singletons['count'] = singletons['count'] / singletons[args.adj_column]

    # convert string genotypes to integers based on config definitino
    replace_dict = config_dict['genotypes']
    geno_asint = geno.replace(replace_dict).replace({1: np.nan})

    # calculate allele frequencies at each site
    ac = np.nansum(geno_asint[samples], axis=1)
    an = np.sum(~np.isnan(geno_asint[samples]), axis=1) * 2
    afs = ac / an

    # only consider sites where allele frequency is between thresholds
    idxs2keep = np.where((afs > 0.1) & (afs < 0.9))[0]
    print ("Using {} genotypes that meet filtering criteria.".format(idxs2keep.shape[0]))
    geno_filtered = geno_asint.iloc[idxs2keep][samples].values
    markers_filtered = geno_asint.iloc[idxs2keep]['marker'].values

    #kinship_matrix = compute_kinship_matrix(geno_filtered)

    # compute the maximum cosine distance between groups of
    # haplotypes at each site in the genotype matrix
    focal_dists = perform_ihd_scan(
        spectra,
        geno_filtered,
    )

    res_df = pd.DataFrame({
        'marker': markers_filtered,
        'distance': focal_dists,
    })
    res_df['k'] = args.k

    # then do permutations
    max_distances = perform_permutation_test(
        spectra,
        geno_filtered,
        n_permutations=args.permutations,
    )

    f, ax = plt.subplots()
    ax.hist(max_distances, bins=20, ec='k', lw=1)

    # compute the 95th percentile of the maximum distance
    # distribution to figure out the distance corresponding to an alpha
    # of 0.05
    for pctile, label in zip((63, 5), ('suggestive', 'significant')):
        score_pctile = np.percentile(max_distances, 100 - pctile)
        res_df[f'{label}_percentile'] = score_pctile
        ax.axvline(x=score_pctile, label=label)

    f.savefig('max_dist.png', dpi=300)

    # compute the median of the maximum distance distribution so
    # that we can compute approximate "odds ratios" for our
    # observed distances
    res_df['null_mean'] = np.mean(max_distances)

    # compute genome-scan adjusted p-values at each marker
    # https://rqtl.org/book/rqtlbook_ch04.pdf, pp113-14
    res_df['pval'] = res_df['distance'].apply(lambda d: np.sum([p >= d for p in max_distances]) / args.permutations)
    res_df.to_csv(args.out, index=False)

if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--singletons")
    p.add_argument("--config")
    p.add_argument("--out")
    p.add_argument("-k", type=int, default=1)
    p.add_argument("-permutations", type=int, default=1_000)
    p.add_argument("-adj_column", default=None)
    p.add_argument("-chrom", default=None)
    args = p.parse_args()

    #numba.set_num_threads(4)


    # lp = LineProfiler()
    # lp.add_function(compute_haplotype_distance)
    # lp_wrapper = lp(main)
    # lp_wrapper(args)
    # lp.print_stats()
    main(args)
