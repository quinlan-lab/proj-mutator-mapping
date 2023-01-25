import pandas as pd
import numpy as np
from typing import List, Tuple
import numba

@numba.njit
def manual_cosine_distance(
    a: np.ndarray,
    b: np.ndarray,
) -> np.float64:
    """function that computes cosine distance between
    two 1d arrays. much faster, since it can be jitted. 

    Args:
        a (np.ndarray): _description_
        b (np.ndarray): _description_

    Returns:
        np.float64: _description_
    """
    dot = a.dot(b)
    a_sumsq, b_sumsq = np.sum(np.square(a)), np.sum(np.square(b))
    a_norm, b_norm = np.sqrt(a_sumsq), np.sqrt(b_sumsq)
    cossim = dot / (a_norm * b_norm)
    return 1 - cossim

@numba.njit
def compute_mean(a: np.ndarray):
    """homebrew function to compute the mean
    on a per-column basis. need to piece this
    out into its own function so that it can be 
    jitted, as numba does not support axis kwargs
    in np.mean

    Args:
        a (np.ndarray): _description_

    Returns:
        _type_: _description_
    """
    empty_a = np.zeros(a.shape[1])
    for i in range(a.shape[1]):
        empty_a[i] = np.mean(a[:, i])
    return empty_a

@numba.njit
def compute_haplotype_distance(
    a_haps: np.ndarray,
    b_haps: np.ndarray,
) -> np.float64:
    """compute the cosine distance between
    two 1d numpy arrays. each array should summarize
    the mutation spectrum in a single group of haplotypes,
    such that the arrays are of shape (M, 1), where M
    is the number of mutations being used to define the spectrum.

    Args:
        a_haps (np.ndarray): _description_
        b_haps (np.ndarray): _description_

    Returns:
        np.float64: cosine distance measure
    """
    # first, sum the spectrum arrays such that we add up the 
    # total number of C>T, C>A, etc. across all samples.
    a_hap_sums = np.sum(a_haps, axis=0)
    b_hap_sums = np.sum(b_haps, axis=0)

    dist = manual_cosine_distance(a_hap_sums, b_hap_sums)

    return dist

@numba.njit
def shuffle_spectra(spectra: np.ndarray, seed: int = 42) -> np.ndarray:
    """shuffle a 2d array of mutation spectrum data so 
    that individual spectra no longer correspond to the
    appropriate sample indices into the rows of that array

    Args:
        spectra (np.ndarray): 2d array of mutation spectrum data
        of shape (S, M), where S is the number of samples and M is the
        number of mutations
        seed (int): random seed

    Returns:
        np.ndarray: shuffled version of the input array
    """
    rng = np.random.default_rng(seed = seed)
    
    idxs = np.arange(spectra.shape[0])
    # shuffle the spectra so that sample idxs no longer
    # correspond to the appropriate spectrum
    rng.shuffle(idxs)
    shuffled_spectra = spectra[idxs, :]
    return shuffled_spectra

@numba.njit()
def permutation_test(
    spectra: np.ndarray,
    genotype_matrix: np.ndarray,
    n_permutations: int = 1_000,
) -> List[np.float64]:
    """conduct a permutation test to assess
    significance of observed psuedo-qtl.

    Args:
        spectra (np.ndarray): 2d array of mutation spectrum data
        of shape (S, M), where S is the number of samples and M is the
        number of mutations
        genotype_matrix (np.ndarray): 2d array of shape (N, S), where
        N is the number of genotyped sites and S is the number of samples.
        n_permutations (int, optional): number of permutations to perform
        (i.e., number of times to shuffle the spectra data and compute cosine
        distances at each marker). Defaults to 1_000.

    Returns:
        List[np.float64]: list of length `n_permutations`, containing the
        maximum cosine distance encountered in each permutation.
    """

    # store max cosdist encountered in each permutation
    max_scores: List[np.float64] = []
    for pi in range(n_permutations):
        if pi > 0 and pi % 100 == 0: print (pi)
        # shuffle the mutation spectra by row
        shuffled_spectra = shuffle_spectra(spectra)

        max_score = 0
        # loop over every site/marker in the genotype matrix
        for ni in range(genotype_matrix.shape[0]):

            # then get the indices of the samples with either
            # genotype at the site. NOTE: we assume that genotypes
            # are coded as 0 or 1.
            a_hap_idxs = np.where(genotype_matrix[ni] == 0)[0]
            b_hap_idxs = np.where(genotype_matrix[ni] == 1)[0]

            a_spectra, b_spectra = (
                shuffled_spectra[a_hap_idxs],
                shuffled_spectra[b_hap_idxs],
            )

            null_dist = compute_haplotype_distance(a_spectra, b_spectra)
            
            # if the distance between haplotype spectra is greater than
            # any we've seen so far in this permutation, store the value
            if null_dist > max_score: max_score = null_dist

        max_scores.append(max_score)

    return max_scores

def compute_spectra(
    singletons: pd.DataFrame,
    k: int = 1,
) -> Tuple[List[str], List[str], np.ndarray,]:
    """compute the mutation spectrum of every sample in 
    the input dataframe of mutation data. we assume the input
    dataframe either contains a single entry for every mutation observed
    in every sample, or contains the aggregate count of each mutation type
    observed in each sample. the input dataframe must have at least three columns:
    one called 'Strain' (denoting the sample), one called 'kmer' (denoting)
    the 3-mer context of the mutation), and 'count' (denoting the number of times
    a mutation of type `kmer` was observed in `Strain`). 

    Args:
        singletons (pd.DataFrame): dataframe containing information about 
        the mutations observed in each sample. the dataframe must have at 
        least three columns: one called 'Strain' (denoting the sample), one 
        called 'kmer' (denoting the 3-mer context of the mutation), and one 
        called 'count' (denoting the number of times a mutation of type `kmer` 
        was observed in `Strain`). 
        k (int, optional): k-mer size of mutations (e.g., k=3 means that 
        we will treat mutations as NXN>NYN, where Ns represent the nucleotide contexts
        on either side of the mutation, whereas k=1 means we will treat
        them as X>Y). Defaults to 1.

    Returns:
        Tuple[List[str], List[str], np.ndarray,]: list of samples in the dataframe 
        (which are indexed in the same order as the rows of the spectra 2d array), 
        list of mutations observed, and a 2d array of mutation spectrum.
    """

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