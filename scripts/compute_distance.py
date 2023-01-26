import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import argparse
import json
import sys
from utils import (
    compute_spectra,
    compute_haplotype_distance,
    permutation_test,
    compute_kinship_diff,
    compute_kinship_matrix,
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

    samples = singletons['Strain'].unique()
    # get the overlap between those and the sample names in the genotype data
    samples_overlap = list(set(samples).intersection(set(geno.columns)))

    # remove known outlier
    samples_overlap = [s for s in samples_overlap if s != "BXD68"]

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

    print (f"Using {len(samples)} samples and {int(np.sum(spectra))} total mutations.")

    # weight spectra if desired
    weight_dict = None
    if args.adj_column:
        # figure out the samples with highest and lowest values of the column
        sorted_by_col = singletons.sort_values(args.adj_column).drop_duplicates("Strain")
        # get each sample's percentage of the max value
        max_colval = np.max(sorted_by_col[args.adj_column].values)
        sorted_by_col['adj_pct'] = sorted_by_col[args.adj_column] / max_colval
        weight_dict = dict(zip(sorted_by_col['Strain'], sorted_by_col['adj_pct']))

    # reweight spectra if desired
    if weight_dict is not None:
        for sample, sample_i in smp2idx.items():
            sample_weight = weight_dict[sample]
            spectra[sample_i] *= sample_weight

    replace_dict = config_dict['genotypes']
    geno_asint = geno.replace(replace_dict).replace({1: np.nan})
    
    # calculate allele frequencies at each site
    ac = np.nansum(geno_asint[samples], axis=1)
    an = np.sum(~np.isnan(geno_asint[samples]), axis=1) * 2
    afs = ac / an

    # only consider sites where allele frequency is between thresholds
    idxs2keep = np.where((afs > 0.05) & (afs < 0.95))[0]
    print ("Using {} genotypes that meet filtering criteria.".format(idxs2keep.shape[0]))
    geno_filtered = geno_asint.iloc[idxs2keep][samples].values
    markers_filtered = geno_asint.iloc[idxs2keep]['marker'].values    

    res = []

    kinship_matrix = compute_kinship_matrix(geno_filtered)

    # loop over every site in the genotype matrix
    for ni in np.arange(geno_filtered.shape[0]):

        # then get the indices of the samples with either
        # genotype at the site
        a_hap_idxs = np.where(geno_filtered[ni] == 0)[0]
        b_hap_idxs = np.where(geno_filtered[ni] == 2)[0]

        a_spectra, b_spectra = (
            spectra[a_hap_idxs],
            spectra[b_hap_idxs],
        )

        focal_dist = compute_haplotype_distance(a_spectra, b_spectra)

        sim_diff = compute_kinship_diff(kinship_matrix, a_hap_idxs, b_hap_idxs)

        res.append({
            'marker': markers_filtered[ni],
            'k': args.k,
            'distance': focal_dist,
            'genetic_difference': sim_diff,
        })

    # then do permutations
    max_distances = permutation_test(
        spectra,
        geno_filtered,
        kinship_matrix,
        n_permutations=args.permutations,
    )

    f, ax = plt.subplots()
    ax.hist(max_distances, bins=20, ec='k', lw=1)

    # compute the 95th percentile of the maximum distance
    # distribution to figure out the distance corresponding to an alpha
    # of 0.05
    res_df = pd.DataFrame(res)
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
    args = p.parse_args()

    #numba.set_num_threads(4)


    # lp = LineProfiler()
    # lp.add_function(compute_haplotype_distance)
    # lp_wrapper = lp(main)
    # lp_wrapper(args)
    # lp.print_stats()
    main(args)
