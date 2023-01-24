import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.metrics.pairwise import cosine_distances, cosine_similarity
from typing import List, Any
import tqdm
import argparse
from line_profiler import LineProfiler
import numba
import re
from collections import Counter
import json 
import sys

@numba.njit()
def permutation_test(
    spectra: np.ndarray,
    genotype_matrix: np.ndarray,
    n_permutations: int = 1_000,
) -> List[np.float64]:

    # then do permutations
    max_scores: List[np.float64] = [] # store max cosdists encountered in each permutation
    for pi in range(n_permutations):
        if pi > 0 and pi % 100 == 0: print (pi)
        idxs = np.arange(spectra.shape[0])
        # shuffle the spectra so that sample idxs no longer
        # correspond to the appropriate spectrum
        np.random.shuffle(idxs)
        shuffled_spectra = spectra[idxs, :]
        max_score = 0
        # loop over every site in the genotype matrix
        for ni in range(genotype_matrix.shape[0]):

            # then get the indices of the samples with either
            # genotype at the site
            a_hap_idxs = np.where(genotype_matrix[ni] == 0)[0]
            b_hap_idxs = np.where(genotype_matrix[ni] == 1)[0]

            a_spectra, b_spectra = (
                shuffled_spectra[a_hap_idxs],
                shuffled_spectra[b_hap_idxs],
            )

            null_dist = compute_haplotype_distance(a_spectra, b_spectra)
            # if null_dist > 0.01:
            #     print (ni)
            if null_dist > max_score: max_score = null_dist

        max_scores.append(max_score)

    return max_scores

@numba.njit
def manual_cosine_distance(a: np.ndarray, b: np.ndarray) -> np.float64:
    dot = a.dot(b)
    a_sumsq, b_sumsq = np.sum(np.square(a)), np.sum(np.square(b))
    a_norm, b_norm = np.sqrt(a_sumsq), np.sqrt(b_sumsq)
    cossim = dot / (a_norm * b_norm)
    return 1 - cossim

@numba.njit
def compute_mean(a: np.ndarray):
    empty_a = np.zeros(a.shape[1])
    for i in range(a.shape[1]):
        empty_a[i] = np.mean(a[:, i])
    return empty_a

@numba.njit
def compute_haplotype_distance(
    a_haps: np.ndarray,
    b_haps: np.ndarray,
) -> np.float64:
    a_hap_sums = np.sum(a_haps, axis=0)
    b_hap_sums = np.sum(b_haps, axis=0)

    dist = manual_cosine_distance(a_hap_sums, b_hap_sums)

    return dist

def generate_null_dist(hap_a, hap_b, n_trials: int = 1_000):
    # concatenate haplotypes
    haps = np.concatenate((hap_a, hap_b))

    dists = []
    for _ in range(n_trials):
        random_mut_hap_idxs = np.random.randint(
            low=0,
            high=haps.shape[0],
            size=hap_a.shape[0],
        )
        random_wt_hap_idxs = np.random.randint(
            low=0,
            high=haps.shape[0],
            size=hap_b.shape[0],
        )

        random_mut_haps = haps[random_mut_hap_idxs]
        random_wt_haps = haps[random_wt_hap_idxs]

        dist = compute_haplotype_distance(
            random_mut_haps,
            random_wt_haps,
        )

        dists.append(dist)

    return dists

def compute_spectra(singletons: pd.DataFrame, k: int = 1):

    # compute 3-mer spectra
    hap_spectra_agg = singletons.groupby(['Strain', 'kmer']).agg({'count': sum}).reset_index()#.rename(columns={0: 'count'})
    # if 1-mer spectra are desired, compute that
    if k == 1:
        # add base mutation type
        hap_spectra_agg['base_mut'] = hap_spectra_agg['kmer'].apply(lambda k: ">".join([k[1], k[5]]))
        hap_spectra_agg = hap_spectra_agg.groupby(['Strain', 'base_mut']).agg({'count': sum}).reset_index()
    # get spectra as per-haplotype vectors of mutation counts
    mut_col = "base_mut" if k == 1 else "kmer"
    spectra = hap_spectra_agg.pivot(index="Strain", columns=mut_col).reset_index().fillna(value=0)
    samples, mutations, spectra = spectra['Strain'].to_list(), [el[1] for el in spectra.columns[1:]], spectra.values[:, 1:]

    return samples, mutations, spectra.astype(np.float32)

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
    geno_asint = geno.replace(replace_dict).replace({-1: np.nan})

    # calculate allele frequencies at each site
    ac = np.nansum(geno_asint[samples], axis=1)
    an = np.sum(~np.isnan(geno_asint[samples]), axis=1)
    afs = ac / an

    # only consider sites where allele frequency is between thresholds
    idxs2keep = np.where((afs > 0.05) & (afs < 0.95))[0]
    print ("Using {} genotypes that meet filtering criteria.".format(idxs2keep.shape[0]))
    geno_filtered = geno_asint.iloc[idxs2keep][samples].values
    markers_filtered = geno_asint.iloc[idxs2keep]['marker'].values

    res = []

    # loop over every site in the genotype matrix
    for ni in np.arange(geno_filtered.shape[0]):

        # then get the indices of the samples with either
        # genotype at the site
        a_hap_idxs = np.where(geno_filtered[ni] == 0)[0]
        b_hap_idxs = np.where(geno_filtered[ni] == 1)[0]

        #print (ni, geno_filtered[ni], np.nansum(geno_filtered[ni]), np.sum(~np.isnan(geno_filtered[ni])))
        a_spectra, b_spectra = (
            spectra[a_hap_idxs],
            spectra[b_hap_idxs],
        )

        focal_dist = compute_haplotype_distance(a_spectra, b_spectra)

        res.append({
            'marker': markers_filtered[ni],
            'k': args.k,
            'distance': focal_dist,
        })

    # then do permutations
    max_scores = permutation_test(
        spectra,
        geno_filtered,
        n_permutations=args.permutations,
    )

    f, ax = plt.subplots()
    ax.hist(max_scores, bins=20, ec='k', lw=1)
    

    # genome-wide significant and suggestive alphas
    res_df = pd.DataFrame(res)
    for pctile, label in zip((63, 5), ('suggestive', 'significant')):
        score_pctile = np.percentile(max_scores, 100 - pctile)
        res_df[f'{label}_percentile'] = score_pctile
        ax.axvline(x=score_pctile, label=label)

    f.savefig('max_dist.png', dpi=300)

    # compute genome-scan adjusted p-values at each marker
    # https://rqtl.org/book/rqtlbook_ch04.pdf, pp113-14
    res_df['pval'] = res_df['distance'].apply(lambda d: np.sum([p >= d for p in max_scores]) / args.permutations)
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
