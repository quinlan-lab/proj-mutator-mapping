import numpy as np
import tqdm
import itertools
import pandas as pd
import numba
from utils import manual_cosine_distance, shuffle_spectra
from typing import List
import copy

class Haplotypes(object):

    def __init__(
        self,
        *,
        n_haplotypes: int,
        frac: float,
        augment: float,
        count: int,
        idx: np.ndarray,
        lambdas: np.ndarray,
    ):
        """initialize the Haplotypes class

        Args:
            n_haplotypes (int): Number of haplotypes to simulate.
            frac (float): Fraction of haplotypes that will contain the mutator.
            augment (float): Effect size of the mutator. The lambda value for \
                the mutation type specified by `idx` will be multiplied by this value.
            count (int): Number of mutations to simulate on each haplotype.
            idx (np.ndarray): Index of the mutation type(s) whose lambda values you \
                 wis to augment.
            lambdas (np.ndarray): Array of lambda values for each mutation type. Presently \
                the array should sum to 1, i.e., the lambda values should correspond to the \
                fraction of de novo mutations expected to belong to each mutation type.
        """
        self.n_haplotypes = n_haplotypes
        self.frac = frac
        self.count = count
        self.augment = augment
        self.idx = idx
        self.lambdas = lambdas

        # figure out the number of wild-type and mutant
        # haplotypes
        self.n_wt_haps = int(self.n_haplotypes * (1 - self.frac))
        self.n_mut_haps = self.n_haplotypes - self.n_wt_haps

    def generate_haplotypes(self, pct_to_augment: float = 1.):
        """create an array of haplotypes of size (N, M), where
        N is the number of haplotypes and M is the number of 
        mutation types.

        Args:
            pct_to_augment (float, optional): The fraction of de novo mutations \
                that should be subject to the effects of the mutator in the mutant \
                haplotypes. If less than 1, the mutant haplotypes will first be assigned \
                (1 - pct_to_augment) of their mutations by taking Poisson draws \
                from the lambda values specified when the class is initialized. Then, \
                those labmdas will be augmented and the remaining `pct_to_augment` mutations \
                will be assigned using the updated lambdas. Defaults to 1..
        """

        # generate "wild-type" haplotypes by taking Poisson
        # draws from the baseline lambda values in each haplotype
        hap_wt = np.random.poisson(lam=np.broadcast_to(
            self.lambdas * self.count, (
                self.n_wt_haps,
                self.lambdas.shape[0],
            )), )

        # generate "mutant" haplotypes by taking Poisson draws from
        # an "augmented" array of lambdas, in which the specified mutation
        # type(s) have been increased.

        # if `pct_to_augment` is not 1 (i.e., if we expect only a fraction
        # of the mutations on haplotypes to be subject to the effects of the
        # mutator), we'll first take Poisson draws from the "wild-type" array of lambdas
        background_count = int(self.count * (1 - pct_to_augment))

        hap_mut = np.random.poisson(lam=np.broadcast_to(
            self.lambdas * background_count, (
                self.n_mut_haps,
                self.lambdas.shape[0],
            )), )

        # then, generate the "focal" mutations under the
        # influence of the mutator
        mut_lambdas = copy.deepcopy(self.lambdas)
        mut_lambdas[self.idx] *= self.augment

        focal_count = self.count - background_count
        hap_mut += np.random.poisson(lam=np.broadcast_to(
            mut_lambdas * focal_count, (
                self.n_mut_haps,
                self.lambdas.shape[0],
            )), )

        self.hap_mut = hap_mut
        self.hap_wt = hap_wt

@numba.njit
def run_permutations(
    *,
    haps: np.ndarray,
    n_markers: int,
    n_permutations: int,
    frac: float,
) -> List[np.float64]:
    """given a simulated set of haplotypes, perform a permutation
    test as follows. in each permutation trial, begin by  shuffling 
    the input haplotypes so that indices into the haplotype array no 
    longer correspond to the correct spectra. then, iterate over a series 
    of toy "markers." at each marker, randomly divide the haplotypes into
    those with the wild-type allele those with the mutant allele. then, compute
    the cosine distance between the aggregate spectra of those two groups. in
    each trial, store the maximum distance encountered at any marker.

    Args:
        haps (np.ndarray): 2D numpy array of haplotypes of shape (N, M), where \
            N is the number of haplotypes and M is the number of mutation types.
        n_markers (int): Number of markers to simulate in each permutation trial.
        n_permutations (int): Number of total permutation to perform.
        frac (float): Expected fraction of haplotypes that should have the WT allele \
            at each marker.

    Returns:
        maximum_distances (List[np.float64]): List of length `n_permutations`, \
            containing the maximum cosine distance encountered at any marker in \
            every permutation trial.
    """
    all_dists = []

    for _ in range(n_permutations):
        # shuffle the haplotypes
        shuffled_spectra = shuffle_spectra(haps)

        max_dist = 0
        # loop over markers
        for _ in range(n_markers):
            # at each marker, choose a random number of samples to make
            # the wt and mut haps, normally distributed around 0.5
            n_mut_haps_to_sample = int(haps.shape[0] * np.random.normal(0.5, 0.05))
            # ensure that we have at least one mutant haplotype
            while n_mut_haps_to_sample < 1:
                n_mut_haps_to_sample = int(haps.shape[0] * np.random.normal(0.5, 0.05))

            n_wt_haps_to_sample = haps.shape[0] - n_mut_haps_to_sample

            # generate indices into the haplotype array. since
            # the spectra are shuffled in each permutation, these indices
            # no longer match the appropriate spectra.
            random_mut_hap_idxs = np.arange(n_mut_haps_to_sample)
            random_wt_hap_idxs = np.arange(n_wt_haps_to_sample) + 1
            random_wt_hap_idxs += np.max(random_mut_hap_idxs)

            # get the aggregate spectra of haplotypes with either allele
            mut_hap_sums = np.sum(shuffled_spectra[random_mut_hap_idxs], axis=0)
            wt_hap_sums = np.sum(shuffled_spectra[random_wt_hap_idxs], axis=0)

            # compute the cosine distance between aggregate spectra
            dist = manual_cosine_distance(mut_hap_sums, wt_hap_sums)
            if dist > max_dist: max_dist = dist

        all_dists.append(max_dist)

    return all_dists


def main():
    rng = np.random.default_rng(42)

    # define the mutation types and expected lambdas to simulate
    base_mutations = ["C>T", "C>A", "C>G", "A>T", "A>C", "A>G"]
    base_lambdas = [0.4, 0.1, 0.075, 0.075, 0.075, 0.275]
    nucs = ["A", "T", "C", "G"]

    kmer_size = 1

    if kmer_size == 1:
        mutations, lambdas = (
            copy.deepcopy(base_mutations),
            copy.deepcopy(base_lambdas),
        )
    else:
        # if we want to do the 3-mer spectrum, we'll just parcel out
        # mutation probabilities equally to every 3-mer associated with
        # a particular "base" mutation type
        mutations, lambdas = [], []
        for m, l in zip(base_mutations, base_lambdas):
            per_k_l = l / 16
            orig, new = m.split('>')
            for fp in nucs:
                for tp in nucs:
                    kmer = f"{fp}{orig}{tp}>{fp}{new}{tp}"
                    mutations.append(kmer)
                    lambdas.append(per_k_l)

    lambdas = np.array(lambdas)
    idx2mut = dict(zip(range(len(mutations)), mutations))

    res = []

    ### ------
    ### DEFINE PARAMETER SEARCH SPACE
    ### ------

    # number of haplotypes to simulate
    n = [100]
    # fraction of samples to add the mutator allele to
    f = [0.5]
    # amounts to augment the mutation probabilities by
    m = [1.01, 1.05, 1.1, 1.2, 1.3, 1.4, 1.5, 2.0]
    # indices of mutations to augment
    if kmer_size == 3:
        # in the case of the 3-mer spectra, we'll augment
        # groups of mutations rather than individual 3-mer types.
        i = [
            list(range(0, 4)),
            list(range(0, 8)),
            list(range(0, 16)),
            list(range(16, 20)),
            list(range(16, 24)),
            list(range(16, 32)),
        ]
    else:
        i = [0, 1, 2]
    # number of total mutations to simulate on each haplotype
    c = [50, 100, 200]
    # fraction of mutations subject to effects of mutator
    p = [1.]
    # number of markers to simulate in each permutation
    k = [1, 10, 100]

    # for each parameter combination, how many replicates to perform
    n_replicates = 10
    # in each replicate, how many permutations to perform
    n_permutations = 100

    for (
            n_haplotypes,
            frac,
            augment,
            count,
            idx,
            pct_to_augment,
            n_markers,
    ) in tqdm.tqdm(itertools.product(n, f, m, c, i, p, k)):

        for rep in range(n_replicates):
            # generate wt and mutant haplotypes
            haplotypes = Haplotypes(
                n_haplotypes=n_haplotypes,
                frac=frac,
                augment=augment,
                count=count,
                idx=np.array(idx),
                lambdas=lambdas,
            )

            haplotypes.generate_haplotypes(pct_to_augment=pct_to_augment)

            haps = np.concatenate((
                haplotypes.hap_mut,
                haplotypes.hap_wt,
            )).astype(np.float32)

            # figure out the indices into the spectra array that
            # correspond to mutant and wt haplotypes
            n_mut_haps = haplotypes.hap_mut.shape[0]
            n_wt_haps = haps.shape[0] - n_mut_haps
            mut_hap_idxs = np.arange(n_mut_haps)
            wt_hap_idxs = np.arange(n_wt_haps) + 1
            wt_hap_idxs += np.max(mut_hap_idxs)

            # compute the cosine distance between the
            # "true" mutant and wt arrays
            focal_dist = manual_cosine_distance(
                np.sum(haps[mut_hap_idxs], axis=0),
                np.sum(haps[wt_hap_idxs], axis=0))

            # simulate
            all_dists = run_permutations(
                haps=haps,
                n_markers=n_markers,
                n_permutations=n_permutations,
                frac=0.5,
            )

            # figure out the mutation types that were
            # changed in each simulation so that we can
            # record them appropriately
            mutation_type = None
            if kmer_size == 1:
                mutation_type = idx2mut[idx]
            else:
                mutation_types = [idx2mut[i] for i in idx]
                base_muts = [">".join([m[1], m[5]]) for m in mutation_types]
                assert len(set(base_muts)) == 1
                base_mut = base_muts[0]
                pct_augmented = int((len(mutation_types) / 16) * 100)
                mutation_type = f"{pct_augmented}% of {base_mut}"

            # calculate the p-value
            p = sum([d >= focal_dist for d in all_dists]) / n_permutations

            res.append({
                    '# of haplotypes': n_haplotypes,
                    '% with mutator': frac,
                    'Mutator effect size': augment,
                    '# of mutations': count,
                    'Mutation type': mutation_type,
                    '# genotyped markers': n_markers,
                    'pval': p,
                    'replicate': rep,
                })

    res_df = pd.DataFrame(res)

    res_df.to_csv('results.csv', index=False)

if __name__ == "__main__":
    # lp = LineProfiler()
    # lp_wrapper = lp(main)
    # lp_wrapper()
    # lp.print_stats()
    main()