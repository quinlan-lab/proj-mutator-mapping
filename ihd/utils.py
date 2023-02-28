import pandas as pd
import numpy as np
from typing import List, Iterable, Sequence, Tuple
import numba
from sklearn.linear_model import LinearRegression
from scipy.stats import chi2_contingency

def chi2_test(a: np.ndarray, b: np.ndarray) -> np.float64:
    observed = np.vstack((a, b))
    #print (observed)
    stat, p, _, _ = chi2_contingency(observed)
    return p


@numba.njit
def manual_cosine_distance(
    a: np.ndarray,
    b: np.ndarray,
) -> np.float64:
    """Compute the cosine distance between
    two 1D numpy arrays. Although methods to compute
    cosine similarity and distance exist in `scipy` and
    `sklearn`, they are not `numba.njit`'able.

    Args:
        a (np.ndarray): A 1D numpy array of size (N, ).
        b (np.ndarray): A 1D numpy array of size (N, ).

    Returns:
        cosine_distance (np.float64): Cosine distance between a and b.
    """
    dot = a.dot(b)
    a_sumsq, b_sumsq = np.sum(np.square(a)), np.sum(np.square(b))
    a_norm, b_norm = np.sqrt(a_sumsq), np.sqrt(b_sumsq)
    cossim = dot / (a_norm * b_norm)
    return 1 - cossim


@numba.njit
def compute_colmean(a: np.ndarray) -> np.ndarray:
    """Compute the mean of a 2D numpy array
    on a per-column basis. Since `numba` does not
    support kwargs in the `np.mean` function, it's
    necessary to piece this out into its own function
    so that it can be decorated with `numba.njit`.

    Args:
        a (np.ndarray): A 2D numpy array of size (N, M).

    Returns:
        colmeans (np.ndarray): A 1D numpy array of size (M, ) containing \
            column-wise means of the input.
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
    """Compute the cosine distance between the aggregate
    mutation spectrum of two collections of haplotype mutation data.
    The input arrays should both be 2D numpy arrays of size (N, M), 
    with rows and columns corresponding to samples and mutation types, 
    respectively. This method will first aggregate the mutation spectra across 
    samples to create two new 1D arrays, each of size (M, ). Then, it 
    will compute the cosine distance between  those two 1D arrays. 

    Args:
        a_haps (np.ndarray): 2D array of size (N, M) containing the mutation \
            spectrum of each sample, where N is the number of samples and M is \
            the number of mutation types.

        b_haps (np.ndarray): 2D array of size (N, M) containing the mutation \
            spectrum of each sample, where N is the number of samples and M \
            is the number of mutation types.

    Returns:
        distance (np.float64): Cosine distance between the aggregate \
            mutation spectra of the two haplotypes.
    """
    # first, sum the spectrum arrays such that we add up the
    # total number of C>T, C>A, etc. across all samples.
    a_hap_sums = np.sum(a_haps, axis=0)
    b_hap_sums = np.sum(b_haps, axis=0)

    dist = manual_cosine_distance(a_hap_sums, b_hap_sums)

    return dist


@numba.njit
def shuffle_spectra(spectra: np.ndarray) -> np.ndarray:
    """Randomly shuffle the rows of a 2D numpy array of 
    mutation spectrum data of size (N, M), where N is the number
    of samples and M is the number of mutation types. Shuffled array
    is returned such that sample mutation spectra no longer correspond to the
    appropriate sample indices into the rows of the array.

    Args:
        spectra (np.ndarray): 2D numpy array of mutation spectrum data of \
        shape (N, M), where N is the number of samples and M is the number \
        of mutation types.

    Returns:
        shuffled_spectra (np.ndarray): The input array, but with shuffled rows.
    """
    idxs = np.arange(spectra.shape[0])
    # shuffle the spectra so that sample idxs no longer
    # correspond to the appropriate spectrum
    np.random.shuffle(idxs)
    shuffled_spectra = spectra[idxs, :]
    return shuffled_spectra


@numba.njit()
def perform_ihd_scan(
    spectra: np.ndarray,
    genotype_matrix: np.ndarray,
) -> List[np.float64]:
    """Iterate over every genotyped marker in the `genotype_matrix`, 
    divide the haplotypes into two groups based on sample genotypes at the 
    marker, and compute the distance between the aggregate mutation spectra
    of each group. Return a list of cosine distances of length (G, ), where
    G is the number of genotyped sites.

    Args:
        spectra (np.ndarray): A 2D numpy array of mutation spectra in \
            all genotyped samples, of size (N, M) where N is the number of \
            samples and M is the number of mutation types.

        genotype_matrix (np.ndarray): A 2D numpy array of genotypes at \
            every genotyped marker, of size (G, N), where G is the number \
            of genotyped sites and N is the number of samples.

    Returns:
        distances (List[np.float64]): List of inter-haplotype cosine \
        distances at every marker.
    """

    # store distances at each marker
    focal_dist: List[np.float32] = []
    # loop over every site in the genotype matrix
    for ni in np.arange(genotype_matrix.shape[0]):
        a_hap_idxs = np.where(genotype_matrix[ni] == 0)[0]
        b_hap_idxs = np.where(genotype_matrix[ni] == 2)[0]

        a_spectra, b_spectra = (
            spectra[a_hap_idxs],
            spectra[b_hap_idxs],
        )

        cur_dist = compute_haplotype_distance(a_spectra, b_spectra)

        focal_dist.append(cur_dist)

    return focal_dist


@numba.njit()
def perform_permutation_test(
    spectra: np.ndarray,
    genotype_matrix: np.ndarray,
    n_permutations: int = 1_000,
    comparison_wide: bool = False,
) -> np.ndarray:
    """Conduct a permutation test to assess the significance of 
    any observed IHD peaks. In each of the `n_permutations` trials, 
    do the following: 1. create a shuffled version of the input mutation `spectra`, so that
    sample names/indices no longer correspond to the appropriate mutation
    spectrum. 2. run an IHD scan by computing the cosine distance between
    the aggregate mutation spectrum of samples with either genotype
    at every marker in the `genotype_matrix`. 3. store the maximum cosine distance encountered at any marker.
    Then, return a list of the maximum cosine distances encountered in each of
    the trials. Alternatively, if `comparison_wide` is True, return a matrix of
    size (P, G), where P is the number of permutations and G is the number of
    genotyped markers, in which we store the cosine distance value encountered at
    every marker in every permutation trial.

    Args:
        spectra (np.ndarray): A 2D numpy array of mutation spectra in all \
            genotyped samples, of size (N, M) where N is the number of samples \
            and M is the number of mutation types.

        genotype_matrix (np.ndarray): A 2D numpy array of genotypes at every \
            genotyped marker, of size (G, N), where G is the number of genotyped \
            sites and N is the number of samples.

        n_permutations (int, optional): Number of permutations to perform \
            (i.e., number of times to shuffle the spectra and compute IHDs at \
            each marker). Defaults to 1_000.

        comparison_wide (bool, optional): Whether to output null distances \
            for each individual marker, as opposed to a genome-wide maximum.
              

    Returns:
        null_distances (np.ndarray): 2D numpy array of size (P, G) \
            where P is the number of permutations and G is either 1 \
            (if we're computing a genome-wide distance threshold) or \
            the number of genotyped markers (if we're computing thresholds \
            at each individual marker).
    """

    # store max cosdist encountered in each permutation, or if desired,
    # per-marker null distances
    n_markers = genotype_matrix.shape[0]
    null_distances: np.ndarray = np.zeros((
        n_permutations,
        n_markers if comparison_wide else 1,
    ))
    
    for pi in range(n_permutations):
        if pi > 0 and pi % 100 == 0: print(pi)
        # shuffle the mutation spectra by row
        shuffled_spectra = shuffle_spectra(spectra)
        # compute the cosine distances at each marker
        perm_distances = perform_ihd_scan(
            shuffled_spectra,
            genotype_matrix,
        )
        if comparison_wide:
            null_distances[pi] = perm_distances 
        else:
            null_distances[pi] = max(perm_distances)

    return null_distances


def compute_spectra(
    mutations: pd.DataFrame,
    k: int = 1,
):
    """Compute the mutation spectrum of every sample in 
    the input pd.DataFrame of mutation data. The input
    dataframe should either contain a single entry for every mutation observed
    in every sample, or the aggregate count of each mutation type
    observed in each sample. The input dataframe must have at least three columns:
    'Strain' (denoting the sample), 'kmer' (denoting the 3-mer context of the 
    mutation -- e.g, CCT>CAT), and 'count' (denoting the number of times
    a mutation of type 'kmer' was observed in 'sample'). If the dataframe contains
    a single entry for every mutation observed in every sample, then the 'count'
    column should contain a value of 1 in every row.

    Args:
        mutations (pd.DataFrame): Pandas dataframe containing information \
            about the mutations observed in each sample. The dataframe must have \
            at least three columns: 'sample' (denoting the sample), 'kmer' (denoting \
            the 3-mer context of the mutation), and 'count' (denoting the number \
            of times a mutation of type 'kmer' was observed in 'sample').

        k (int, optional): k-mer size of mutations (e.g., k=3 means that we \
            will treat mutations as NXN NYN, where Ns represent the nucleotide \
            contexts on either side of the mutation, whereas k=1 means we will \
            treat them as X Y). Defaults to 1.

    Returns:
        samples (List[str]): A list of samples in the dataframe (which are \
            indexed in the same order as the rows of the 2D `spectra` array).
        mutation_types (List[str]): A list of the unique mutation types \
            observed in the dataframe.
        spectra (np.ndarray): A 2D numpy array of mutation spectra in \
            all genotyped samples, of size (N, M) where N is the number \
            of samples and M is the number of mutation types.
    """

    # compute 3-mer spectra
    hap_spectra_agg = mutations.groupby(['sample', 'kmer']).agg({
        'count': sum
    }).reset_index()
    # if 1-mer spectra are desired, compute that
    if k == 1:
        # add base mutation type
        hap_spectra_agg['base_mut'] = hap_spectra_agg['kmer'].apply(
            lambda k: ">".join([k[1], k[5]]))
        hap_spectra_agg = hap_spectra_agg.groupby(['sample', 'base_mut']).agg({
            'count':
            sum
        }).reset_index()
    # get spectra as per-haplotype vectors of mutation counts
    mut_col = "base_mut" if k == 1 else "kmer"
    spectra = hap_spectra_agg.pivot(
        index="sample", columns=mut_col).reset_index().fillna(value=0)
    samples, mutations, spectra = spectra['sample'].to_list(), [
        el[1] for el in spectra.columns[1:]
    ], spectra.values[:, 1:]

    return samples, mutations, spectra.astype(np.float64)
