import numpy as np
import tqdm
import itertools
import pandas as pd
import numba
from typing import List
from utils import (
    perform_ihd_scan,
    perform_permutation_test,
    compute_genotype_similarity,
)
import argparse

@numba.njit
def simulate_genotypes(
    n_markers: int,
    n_haplotypes: int,
    exp_af: float = 0.5,
) -> np.ndarray:
    """Simulate genotypes on haplotypes. 

    Args:
        n_markers (int): Number of sites at which to simulate genotypes.
        
        n_haplotypes (int): Number of haplotypes to simulate.

        exp_af (float, optional): Expected frequency of ALT allele at each site. Defaults to 0.5.

    Returns:
        np.ndarray: Matrix of size (g, n), where g is the number of markers \
            and n is the number of haplotypes.
    """
    genotype_matrix = np.zeros((n_markers, n_haplotypes))

    # simulate genotypes at a specific number of markers.
    # for each haplotype, take a random draw from a uniform
    # distribution at each marker. if that draw is <= the
    # expected marker allele frequency, assign that marker's
    # genotype to be 1. otherwise, keep it at 0.
    for haplotype_idx in np.arange(n_haplotypes):
        genotype_probs = np.random.random(n_markers)
        # convert probs to genotypes
        alt_alleles = np.where(genotype_probs <= exp_af)[0]
        genotype_matrix[alt_alleles, haplotype_idx] = 2
    return genotype_matrix


@numba.njit
def poisson_from_arr(lambdas: np.ndarray) -> np.ndarray:
    """Take Poisson draws from an array of lambda values
     that can be multi-dimensional. Numba doesn't allow
     lambdas to be arrays in np.random.poisson, so we have
     to wrap it up and njit it.

    Args:
        lambdas (np.ndarray): Any-dimensional array of lambdas.

    Returns:
        np.ndarray: Any-dimensional array of counts, of the same \
            shape as the input lambda array.
    """
    output_arr = np.zeros(lambdas.shape)
    for i in np.arange(lambdas.shape[0]):
        for j in np.arange(lambdas.shape[1]):
            sample = np.random.poisson(lambdas[i, j])
            output_arr[i, j] = sample
    return output_arr

@numba.njit
def run_simulation_trials(
    base_lambdas: np.ndarray,
    genotype_matrix: np.ndarray,
    mutation_type_idx: int = 0,
    effect_size: float = 1.1,
    n_permutations: int = 100,
    n_mutations: int = 50,
    n_haplotypes: int = 100,
    n_markers: int = 1000,
    number_of_trials: int = 100,
    f_with_mutator: float = 1.,
) -> List[np.float64]:
    """Run a selected number of simulation trials.

    Args:
        base_lambdas (np.ndarray): 1D numpy array of lambda values corresponding to the \
            expected number of each k-mer mutation type we observe on a haplotype.

        genotype_matrix (np.ndarray): 2D numpy array of size (g, n), where g is the number \
            of genotyped sites and n is the number of simulated haplotypes.

        mutation_type_idx (int, optional): Index of the mutation type to augment by the specified \
            effect size. Defaults to 0.

        effect_size (float, optional): Effect size of the mutator allele. Defaults to 1.1.

        n_permutations (int, optional): Number of permutations to use when establishing distance \
            thresholds for significance. Defaults to 100.

        n_mutations (int, optional): Total number of mutations to simulate on each haplotype. Defaults to 50.

        n_haplotypes (int, optional): Number of haplotypes to simulate. Defaults to 100.

        n_markers (int, optional): Number of sites at which to simulate genotypes. Defaults to 1000.

        number_of_trials (int, optional): Number of trials to run for every combination of \
            simulation parameters. Defaults to 100.

        f_with_mutator (float, optional): The fraction of haplotypes with the ALT allele at the \
            simulated mutator locus that actually possess the mutator allele and augmented mutation \
            spectra. I.e., the degree to which the ALT allele at the mutator locus truly tags \
            the mutator allele. Defaults to 1..

    Returns:
        List[np.float64]: List of p-values from `number_of_trials` trials.
    """

    pvals = []

    for _ in range(number_of_trials):

        # create a matrix of lambda values
        lambdas = np.broadcast_to(
            base_lambdas * n_mutations,
            (
                n_haplotypes,
                base_lambdas.shape[0],
            ),
        ).copy()

        # pick a marker at which we'll artificially put an
        # "association" between genotypes and mutation spectra.
        focal_marker = np.random.randint(0, n_markers)
        # at that marker, get the haplotype indices with alt alleles
        alt_haplotypes = np.where(genotype_matrix[focal_marker] == 2)[0]

        # figure out how many of the haplotypes at the "mutator locus"
        # actually carry the mutator allele. in other words, what if all of
        # the markers have AF = 0.5, but the "focal" marker *imperfectly tags* the mutator
        # locus? such that 1/2 of haplotypes carry the marker, but only a fraction of those
        # actually carry the true mutator allele.
        n_true_alt_haplotypes = int(alt_haplotypes.shape[0] * f_with_mutator)
        alt_haplotypes = np.random.choice(alt_haplotypes, n_true_alt_haplotypes)

        # ensure that at least one haplotype at this site has the
        # alternate allele. otherwise, return a p-value of 1.
        if alt_haplotypes.shape[0] < 1:
            pvals.append(1.)
            continue

        # augment the lambda on the alt haplotypes by an effect size
        for ai in alt_haplotypes:
            lambdas[ai, mutation_type_idx] *= effect_size

        # simulate mutation spectra using new lambdas
        mutation_spectra = poisson_from_arr(lambdas)

        genotype_similarity = compute_genotype_similarity(genotype_matrix)

        # run an IHD scan
        focal_dists = perform_ihd_scan(
            mutation_spectra,
            genotype_matrix,
            genotype_similarity,
        )

        # and get null
        null_distances = perform_permutation_test(
            mutation_spectra,
            genotype_matrix,
            genotype_similarity,
            n_permutations=n_permutations,
        )

        pval = np.sum(null_distances >= focal_dists[focal_marker]) / n_permutations
        pvals.append(pval)

    return pvals


def main(args):

    # define parameter space for simulations
    number_of_markers = [1_000]
    number_of_haplotypes = [50, 100]
    number_of_mutations = [20, 100, 500]
    number_of_permutations = [100]
    mutation_types = ["C>T", "C>A", "C>G"]
    effect_sizes = list(np.arange(1, 1.5, 0.1))
    expected_marker_afs = [0.5]
    number_of_trials = 100
    tag_strengths = [1.]

    base_mutations = ["C>T", "CpG>TpG", "C>A", "C>G", "A>T", "A>C", "A>G"]
    base_lambdas = np.array([0.29, 0.17, 0.12, 0.075, 0.1, 0.075, 0.17])
    mut2idx = dict(zip(base_mutations, range(len(base_mutations))))

    res_df = []

    for (
            n_markers,
            n_haplotypes,
            n_mutations,
            n_permutations,
            mutation_type,
            effect_size,
            exp_af,
            tag_strength,
    ) in tqdm.tqdm(
            itertools.product(
                number_of_markers,
                number_of_haplotypes,
                number_of_mutations,
                number_of_permutations,
                mutation_types,
                effect_sizes,
                expected_marker_afs,
                tag_strengths,
            )):

        mutation_type_idx = mut2idx[mutation_type]

        # NOTE: simulate genotypes separately in each trial? or just
        # simulate the mutation spectra independently?
        genotype_matrix = simulate_genotypes(n_markers, n_haplotypes, exp_af=exp_af)

        trial_pvals = run_simulation_trials(
            base_lambdas,
            genotype_matrix,
            mutation_type_idx=mutation_type_idx,
            effect_size=effect_size,
            n_permutations=n_permutations,
            n_mutations=n_mutations,
            n_haplotypes=n_haplotypes,
            n_markers=n_markers,
            number_of_trials=number_of_trials,
            f_with_mutator=tag_strength,
        )

        for i, pval in enumerate(trial_pvals):
            res_df.append({
                "n_markers": n_markers,
                "n_haplotypes": n_haplotypes,
                "n_permutations": n_permutations,
                "n_mutations": n_mutations,
                "effect_size": effect_size,
                "mutation_type": mutation_type,
                "exp_af": exp_af,
                "tag_strength": tag_strength,
                "trial": i,
                "pval": pval,
            })

    res_df = pd.DataFrame(res_df)
    res_df.to_csv(args.out, index=False)

if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--out", help="""path to output file""")
    args = p.parse_args()
    main(args)