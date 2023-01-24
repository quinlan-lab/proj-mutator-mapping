import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.metrics.pairwise import cosine_distances
from typing import List
import tqdm
import argparse
from compute_distance import compute_spectra, generate_null_dist


def main(args):

    geno = pd.read_csv(args.geno, sep="\t", skiprows=21)

    metadata = pd.read_excel(args.metadata)[[
        "GeneNetwork name",
        "bam_name",
        "true_epoch",
    ]].dropna()

    gn2bam = dict(zip(metadata['GeneNetwork name'], metadata['bam_name']))

    singletons = pd.read_csv(args.singletons)


    # for the null permutation test, shuffle the rows of the spectra
    # dataframe every time. otherwise keep it the same.
    samples, spectra = compute_spectra(singletons, k=args.k)
    smp2idx = dict(zip(samples, range(len(samples))))

    res = []

    for n_samps in tqdm.tqdm(range(1, len(samples))):

        # randomly grab B and D haplotype samples
        b_hap_samples = np.random.choice(samples, size=n_samps)
        d_hap_samples = [s for s in samples if s not in b_hap_samples]

        b_hap_idxs = [smp2idx[s] for s in b_hap_samples]
        d_hap_idxs = [smp2idx[s] for s in d_hap_samples]

        b_spectra = spectra[b_hap_idxs]
        d_spectra = spectra[d_hap_idxs]

        null_dists = generate_null_dist(
            b_spectra,
            d_spectra,
            n_trials=args.n_trials,
        )

        res.append({
            'n_focal_haps': n_samps,
            'k': args.k,
            'null_median': np.median(null_dists),
        })

    res_df = pd.DataFrame(res)
    res_df.to_csv(args.out, index=False)

if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--geno")
    p.add_argument("--metadata")
    p.add_argument("--singletons")
    p.add_argument("--out")
    p.add_argument("-n_trials", type=int, default=1_000)
    p.add_argument("-k", type=int, default=1)
    args = p.parse_args()
    main(args)
