import numpy as np
import tqdm
import itertools
import pandas as pd
from scipy.spatial.distance import cosine as cosdist
import numba
import scipy.stats as ss
import sys

# adding Folder_2 to the system path
sys.path.insert(0, '/Users/tomsasani/quinlanlab/proj-mutator-mapping/')
from ihd.utils import compute_manual_chisquare
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
        self.n_haplotypes = n_haplotypes
        self.frac = frac
        self.count = count
        self.augment = augment
        self.idx = idx
        self.lambdas = lambdas

        # get the number of WT haplotypes by taking a poisson draw
        # from the total number of haplotypes, using the expected fraction
        # of WT haplotypes as the lambda
        self.n_wt_haps = int(np.random.poisson(
            int(self.n_haplotypes * (1 - self.frac)),
            size=1,
        )[0])
        # and get the number of MUT haplotypes
        self.n_mut_haps = self.n_haplotypes - self.n_wt_haps
        # adjust counts if necessary
        while self.n_mut_haps < 1: self.n_mut_haps += 1
        while self.n_wt_haps < 1: self.n_wt_haps += 1


    def generate_haplotypes(self, pct_to_augment: float = 1.):

        # generate "wild-type" haplotypes
        hap_wt = np.random.poisson(lam=np.broadcast_to(
            self.lambdas * self.count,
            (
                self.n_wt_haps,
                self.lambdas.shape[0],
            ),
        ), )

        # and generate "mutant" haplotypes.
        # if we are augmenting *all* of the mutations in these
        # haplotypes (i.e., if `pct_to_augment` is set to 1), we're
        # assuming that each of these mutations is a DNM
        # that is being affected by the mutator allele.
        # if `pct_to_augment` is not 1, then we're assuming that each
        # haplotype has a background set of mutations that were already
        # present on the haplotype before the mutator exerted its effect.

        background_count = int(self.count * (1 - pct_to_augment))

        hap_mut = np.random.poisson(lam=np.broadcast_to(
            self.lambdas * background_count,
            (
                self.n_mut_haps,
                self.lambdas.shape[0],
            ),
        ), )

        # then, generate the "focal" mutations under the influence of
        # the mutator
        mut_lambdas = self.lambdas.copy()
        mut_lambdas[self.idx] *= self.augment
        focal_count = self.count - background_count
        hap_mut += np.random.poisson(lam=np.broadcast_to(
            mut_lambdas * focal_count,
            (
                self.n_mut_haps,
                self.lambdas.shape[0],
            ),
        ), )
        self.hap_mut = hap_mut
        self.hap_wt = hap_wt


@numba.njit
def run_permutations(
    *,
    wt_haps: np.ndarray,
    mut_haps: np.ndarray,
    n_markers: int,
    n_permutations: int,
    frac: float,
):
    all_dists = []

    haps = np.concatenate((wt_haps, mut_haps))

    for _ in range(n_permutations):
        # shuffle the haplotypes
        shuffled_idxs = np.arange(haps.shape[0])
        np.random.shuffle(shuffled_idxs)
        shuffled_spectra = haps[shuffled_idxs].astype(np.float32)

        max_dist = 0
        n_mut_haps_to_sample = mut_haps.shape[0]
        n_wt_haps_to_sample = haps.shape[0] - n_mut_haps_to_sample

        for _ in range(n_markers):

            # at each marker, grab a random fraction of WT and MUT haplotype indices
            random_mut_hap_idxs = np.arange(n_mut_haps_to_sample)
            random_wt_hap_idxs = np.arange(n_wt_haps_to_sample) + 1
            random_wt_hap_idxs += np.max(random_mut_hap_idxs)

            random_mut_haps = shuffled_spectra[random_mut_hap_idxs]
            random_wt_haps = shuffled_spectra[random_wt_hap_idxs]

            a_hap_sums = np.sum(random_mut_haps, axis=0)
            b_hap_sums = np.sum(random_wt_haps, axis=0)
            
            dist = compute_manual_chisquare(a_hap_sums, b_hap_sums)
            
            if dist > max_dist: max_dist = dist

        all_dists.append(max_dist)

    return all_dists


def main():

    base_mutations = ["C>T", "CpG>TpG", "C>A", "C>G", "A>T", "A>C", "A>G"]
    base_mutations = [m.replace(">", r"$\rightarrow$") for m in base_mutations]
    base_lambdas = [0.29, 0.17, 0.12, 0.075, 0.1, 0.075, 0.17]
    nucs = ["A", "T", "C", "G"]

    kmer_size = 1

    if kmer_size == 1:
        mutations, lambdas = base_mutations.copy(), base_lambdas.copy()
    else:
        # if we want to do the 3-mer spectrum, we'll just parcel out
        # mutation probabilities equally to every 3-mer associated with
        # a particular "base" mutation type
        mutations, lambdas = [], []
        for m, l in zip(base_mutations, base_lambdas):
            if m == "C>T":
                l = 0.46
            if m == "CpG>TpG": continue
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

    n_haplotypes = [50, 100]  # number of haplotypes to simulate
    frac = [0.5]  # fraction of samples to add a mutator to
    effect_size = list(np.arange(1, 1.5, 0.05))  # amount to augment the mutation probability (lambda) by
    if kmer_size == 3:
        mutation_idxs = [
            list(range(0, 4)),
            list(range(0, 8)),
            list(range(0, 16)),
            list(range(16, 20)),
            list(range(16, 24)),
            list(range(16, 32)),
        ]
    else:
        mutation_idxs = [0, 2, 3]
    mutation_count = [10, 50, 100, 500]  # number of mutations to simulate per haplotypes
    pct_to_augment = [1.]  # fraction of mutations subject to effects of mutator
    n_markers = [1]  # number of markers used

    for (
            n,
            f,
            m,
            c,
            i,
            p,
            k,
    ) in tqdm.tqdm(
            itertools.product(
                n_haplotypes,
                frac,
                effect_size,
                mutation_count,
                mutation_idxs,
                pct_to_augment,
                n_markers,
            )):

        replicates = 100
        n_permutations = 1_000

        for rep in range(replicates):
            # generate wt and mutant haplotypes
            haplotypes = Haplotypes(
                n_haplotypes=n,
                frac=f,
                augment=m,
                count=c,
                idx=np.array(i),
                lambdas=lambdas,
            )

            haplotypes.generate_haplotypes(pct_to_augment=p)

            # in the GWAS scenario, we are using markers that are approximately
            # 50% AF. so, when trying to detect our mutator, the best case scenario
            # would be testing at a marker where all of the mutant haplotypes have the
            # ALT genotype and all of the WT haplotypes have the WT genotype (or vice versa).
            # but if the mutator is rare, then a bunch of WT haplotypes will be lumped in and
            # also have the ALT genotype.
            haps = np.concatenate((
                haplotypes.hap_mut,
                haplotypes.hap_wt,
            )).astype(np.float32)

            # generate a 50/50 split of mut and wt haplotypes, with the expectation
            # that all of the mut haplotypes are in the same group at this marker
            n_mut_haps = haplotypes.hap_mut.shape[0]
            n_wt_haps = haps.shape[0] - n_mut_haps
            mut_hap_idxs = np.arange(n_mut_haps)
            wt_hap_idxs = np.arange(n_wt_haps) + 1
            wt_hap_idxs += np.max(mut_hap_idxs)

            focal_dist = compute_manual_chisquare(
                np.sum(haps[mut_hap_idxs], axis=0),
                np.sum(haps[wt_hap_idxs], axis=0))


            # simulate with expectation that all markers are at 50%
            all_dists = run_permutations(
                wt_haps=haplotypes.hap_wt,
                mut_haps=haplotypes.hap_mut,
                n_markers=k,
                n_permutations=n_permutations,
                frac=0.5,
            )

            mutation_type = None
            if kmer_size == 1:
                mutation_type = idx2mut[i]
            else:
                mutation_types = [idx2mut[idx] for idx in i]
                base_muts = [">".join([m[1], m[5]]) for m in mutation_types]
                assert len(set(base_muts)) == 1
                base_mut = base_muts[0]
                pct_augmented = int((len(mutation_types) / 16) * 100)
                mutation_type = f"{pct_augmented}% of {base_mut}"

            res.append({
                    '# of haplotypes': n,
                    '% with mutator': f,
                    'Mutator effect size': m,
                    '# of mutations': c,
                    'Mutation type': mutation_type,
                    '# genotyped markers': k,
                    'pval': sum([d >= focal_dist for d in all_dists]) / n_permutations,
                    'replicate': rep,
                })

    res_df = pd.DataFrame(res)

    res_df.to_csv('results.csv', index=False)

if __name__ == "__main__":
    main()