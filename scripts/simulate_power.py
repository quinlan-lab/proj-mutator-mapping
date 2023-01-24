import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.metrics.pairwise import cosine_distances, cosine_similarity
import tqdm
import itertools
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from sklearn.cluster import AgglomerativeClustering, KMeans
import hdbscan
from sklearn.metrics import silhouette_score
from scipy.spatial.distance import cosine as cosdist

PROJDIR = "/scratch/ucgd/lustre-work/quinlan/u1006375/proj-human-spectra"

class Haplotypes(object):

    def __init__(
        self,
        *,
        n_haplotypes: int,
        frac: float,
        augment: float,
        count: int,
        idx: int,
        lambdas: np.ndarray,
    ):
        self.n_haplotypes = n_haplotypes
        self.frac = frac
        self.count = count
        self.augment = augment
        self.idx = idx
        self.lambdas = lambdas

        # get the number of wt haplotypes by taking a poisson draw
        # from the expected ratio of 50/50
        self.n_wt_haps = int(np.random.poisson(
            int(self.n_haplotypes * (1 - self.frac)),
            size=1,
        )[0])
        self.n_mut_haps = self.n_haplotypes - self.n_wt_haps

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

    def calculate_empirical_p(self, n_trials: int = 1_000):
        self.generate_focal_dist()
        self.generate_null_dist(n_trials = n_trials)

        p = np.sum([d >= self.focal_dist for d in self.null_dists])

        return p / n_trials

rng = np.random.default_rng(42)

mutations = ["C>T", "C>A", "C>G", "A>T", "A>C", "A>G"]
idx2mut = dict(zip(range(len(mutations)), mutations))
lambdas = np.array([0.4, 0.1, 0.075, 0.075, 0.075, 0.275])

res = []

n = [700] # number of haplotypes to simulate
f = [0.5] # fraction of samples to add an eQTL to
m = [1.1, 1.5, 2] # amount to augment the mutation probability (lambda) by
i = list(range(6)) # mutation to augment
c = [8] # number of mutations to simulate per haplotypes
p = [1.] # fraction of mutations subject to effects of mutator
k = [1, 10, 100] # number of markers used

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


    # generate wt and mutant haplotypes
    haplotypes = Haplotypes(
        n_haplotypes=n_haplotypes,
        frac=frac,
        augment=augment,
        count=count,
        idx=idx,
        lambdas=lambdas,
    )

    haplotypes.generate_haplotypes(pct_to_augment=pct_to_augment)
    haplotypes.generate_focal_dist()

    haps = np.concatenate((haplotypes.hap_mut, haplotypes.hap_wt))

    n_permutations = 1_000
    all_dists = []

    for p in range(n_permutations):
        # shuffle the haplotypes 
        shuffled_idxs = np.arange(haps.shape[0])
        np.random.shuffle(shuffled_idxs)
        shuffled_spectra = haps[shuffled_idxs]

        max_dist = 0

        for marker in range(n_markers):

            n_haps_to_sample = np.random.poisson(int(haps.shape[0] * frac))

            # at each marker, grab a random fraction of WT and MUT haplotype indices
            random_mut_hap_idxs = np.random.randint(
                            low=0,
                            high=haps.shape[0],
                            size=n_haps_to_sample,
                        )
            random_wt_hap_idxs = shuffled_idxs[~np.in1d(shuffled_idxs, random_mut_hap_idxs)]

            random_mut_haps = haps[random_mut_hap_idxs]
            random_wt_haps = haps[random_wt_hap_idxs]

            a_hap_sums = np.sum(random_mut_haps, axis=0)
            b_hap_sums = np.sum(random_wt_haps, axis=0)

            dist = cosdist(a_hap_sums, b_hap_sums)

            if dist > max_dist: max_dist = dist

        all_dists.append(max_dist)

    # what fraction of randomly shuffled "marker haplotypes" produce
    # a cosine distance >= the focal distance
    res.append({
            'n_haplotypes': n_haplotypes,
            'frac_with_eQTL': frac,
            'augment_factor': augment,
            'mutation_count': count,
            'mutation': idx2mut[idx],
            'pct_to_augment': pct_to_augment,
            'n_markers': n_markers,
            'pval': sum([d >= haplotypes.focal_dist for d in all_dists]) / n_permutations,
        })

res_df = pd.DataFrame(res)


res_df.to_csv('results.csv', index=False)