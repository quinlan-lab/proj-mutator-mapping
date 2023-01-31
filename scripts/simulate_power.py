import numpy as np
import tqdm
import itertools
import pandas as pd
from scipy.spatial.distance import cosine as cosdist
import numba

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

        # get the number of wt haplotypes by taking a poisson draw
        # from the expected ratio
        self.n_wt_haps = int(np.random.poisson(
            int(self.n_haplotypes * (1 - self.frac)),
            size=1,
        )[0])
        self.n_mut_haps = int(np.random.poisson(
            int(self.n_haplotypes * self.frac),
            size=1,
        )[0])

        if self.n_mut_haps < 1: self.n_mut_haps += 1
        if self.n_wt_haps < 1: self.n_wt_haps += 1

    def compute_haplotype_distance(self, a_haps: np.ndarray, b_haps: np.ndarray,):
        a_hap_sums = np.sum(a_haps, axis=0)
        b_hap_sums = np.sum(b_haps, axis=0)

        dist = cosdist(a_hap_sums, b_hap_sums)

        return dist

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
        # assuming that each of these mutations is effectively a DNM
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

    def generate_focal_dist(self):
        self.focal_dist = self.compute_haplotype_distance(
            self.hap_mut,
            self.hap_wt,
        )

@numba.njit
def manual_cosine_distance(a: np.ndarray, b: np.ndarray) -> np.float64:
    dot = a.dot(b)
    a_sumsq, b_sumsq = np.sum(np.square(a)), np.sum(np.square(b))
    a_norm, b_norm = np.sqrt(a_sumsq), np.sqrt(b_sumsq)
    cossim = dot / (a_norm * b_norm)
    return 1 - cossim

@numba.njit
def run_permutations(
    *,
    haps: np.ndarray,
    n_markers: int,
    n_permutations: int,
    frac: float,
):
    all_dists = []

    for _ in range(n_permutations):
        # shuffle the haplotypes
        shuffled_idxs = np.arange(haps.shape[0])
        np.random.shuffle(shuffled_idxs)
        shuffled_spectra = haps[shuffled_idxs].astype(np.float32)

        max_dist = 0

        for _ in range(n_markers):

            # at each marker, choose a random number of samples to make
            # the wt and mut haps, poisson distributed around 0.5
            n_mut_haps_to_sample = np.random.poisson(int(haps.shape[0] * frac))
            # ensure that we're relative close to a 50/50 balance
            while (n_mut_haps_to_sample / haps.shape[0]) < 0.25 or (n_mut_haps_to_sample / haps.shape[0]) > 0.75:
                n_mut_haps_to_sample = np.random.poisson(int(haps.shape[0] * frac))

            n_wt_haps_to_sample = haps.shape[0] - n_mut_haps_to_sample
            # at each marker, grab a random fraction of WT and MUT haplotype indices
            random_mut_hap_idxs = np.arange(n_mut_haps_to_sample)
            random_wt_hap_idxs = np.arange(n_wt_haps_to_sample) + 1
            random_wt_hap_idxs += np.max(random_mut_hap_idxs)

            random_mut_haps = shuffled_spectra[random_mut_hap_idxs]
            random_wt_haps = shuffled_spectra[random_wt_hap_idxs]

            a_hap_sums = np.sum(random_mut_haps, axis=0)
            b_hap_sums = np.sum(random_wt_haps, axis=0)

            dist = manual_cosine_distance(a_hap_sums, b_hap_sums)

            if dist > max_dist: max_dist = dist

        all_dists.append(max_dist)

    return all_dists


def main():
    rng = np.random.default_rng(42)

    base_mutations = ["C>T", "C>A", "C>G", "A>T", "A>C", "A>G"]
    base_lambdas = [0.4, 0.1, 0.075, 0.075, 0.075, 0.275]
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
            per_k_l = l / 16
            orig, new = m.split('>')
            for fp in nucs:
                for tp in nucs:
                    kmer = f"{fp}{orig}{tp}>{fp}{new}{tp}"
                    mutations.append(kmer)
                    lambdas.append(per_k_l)

    lambdas = np.array(lambdas)
    idx2mut = dict(zip(range(len(mutations)), mutations))
    mut2idx = {v:k for k,v in idx2mut.items()}

    res = []

    n = [50, 100] # number of haplotypes to simulate
    f = [0.5] # fraction of samples to add a mutator to
    m = [1.01, 1.05, 1.1, 1.2, 1.3, 1.4, 1.5, 2.0] # amount to augment the mutation probability (lambda) by
    if kmer_size == 3:
        i = [[21, 24, 28, 29, 44]] # mutation(s) to augment
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
    c = [50, 100, 200] # number of mutations to simulate per haplotypes
    p = [1.] # fraction of mutations subject to effects of mutator
    k = [1] # number of markers used

    for (
            n_haplotypes,
            frac,
            augment,
            count,
            idx,
            pct_to_augment,
            n_markers,
    ) in tqdm.tqdm(itertools.product(
            n,
            f,
            m,
            c,
            i,
            p,
            k,
    )):

        replicates = 100
        n_permutations = 500

        for rep in range(replicates):
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

            # in the GWAS scenario, we are using markers that are approximately
            # 50% AF. so, when trying to detect our mutator, the best case scenario
            # would be testing at a marker where all of the mutant haplotypes have the
            # ALT genotype and all of the WT haplotypes have the WT genotype (or vice versa).
            # but if the mutator is rare, then a bunch of WT haplotypes will be lumped in and
            # also have the ALT genotype.
            haps = np.concatenate((haplotypes.hap_mut, haplotypes.hap_wt)).astype(np.float32)

            # generate a 50/50 split of mut and wt haplotypes, with the expectation
            # that all of the mut haplotypes are in the same group at this marker
            n_mut_haps = haplotypes.hap_mut.shape[0]
            n_wt_haps = haps.shape[0] - n_mut_haps
            mut_hap_idxs = np.arange(n_mut_haps)
            wt_hap_idxs = np.arange(n_wt_haps) + 1
            wt_hap_idxs += np.max(mut_hap_idxs)

            focal_dist = manual_cosine_distance(
                np.sum(haps[mut_hap_idxs], axis=0),
                np.sum(haps[wt_hap_idxs], axis=0))

            # simulate with expectation that all markers are at 50%
            all_dists = run_permutations(
                haps=haps,
                n_markers=n_markers,
                n_permutations=n_permutations,
                frac=0.5,
            )

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

            res.append({
                    '# of haplotypes': n_haplotypes,
                    '% with mutator': frac,
                    'Mutator effect size': augment,
                    '# of mutations': count,
                    'Mutation type': mutation_type,
                    '# genotyped markers': n_markers,
                    'pval': sum([d >= focal_dist for d in all_dists]) / n_permutations,
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