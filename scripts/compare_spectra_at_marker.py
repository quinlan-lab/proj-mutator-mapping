import pandas as pd
import argparse
from compute_distance import compute_spectra
import numpy as np
from compare_mutation_spectra import mutation_comparison
import matplotlib.pyplot as plt 
import seaborn as sns 
from sklearn.metrics.pairwise import cosine_distances, cosine_similarity

PROJDIR = "/scratch/ucgd/lustre-work/quinlan/u1006375/proj-mutator-mapping"

geno = pd.read_csv(f"{PROJDIR}/data/genotypes/bxd.geno").set_index("marker")

marker = "rs46183100" 

chrom = "chr17"

k = 1

geno_at_marker = geno.loc[marker]

singletons = pd.read_csv(f"{PROJDIR}/data/singletons/bxd/annotated_filtered_singletons.csv")
singletons['count'] = 1
#singletons = singletons[singletons['chrom'] == chrom]

samples, mutations, spectra = compute_spectra(singletons, k=k)
smp2idx = dict(zip(samples, range(len(samples))))
mut2idx = dict(zip(mutations, range(len(mutations))))
samples_shared = list(set(samples).intersection(set(geno.columns[1:])))

a_smp_i = np.array([smp2idx[k] for k,v in geno_at_marker.to_dict().items() if k in samples_shared and v == "B"])
b_smp_i = np.array([smp2idx[k] for k,v in geno_at_marker.to_dict().items() if k in samples_shared and v == "D"])

a_spectra, b_spectra = spectra[a_smp_i], spectra[b_smp_i]
all_spectra = np.concatenate((a_spectra, b_spectra))
all_spectra_fracs = all_spectra / np.sum(all_spectra, axis=1)[:, np.newaxis]

a_spectra_sum = np.sum(spectra[a_smp_i], axis=0)
b_spectra_sum = np.sum(spectra[b_smp_i], axis=0)

if k == 3:
    mutation_comparison(b_spectra_sum, a_spectra_sum, mut2idx)
elif k == 1:
    f, ax = plt.subplots()
    spectra_fracs = spectra / np.sum(spectra, axis=1)[:, np.newaxis]
    df = []
    for si, s in enumerate(samples):
        for m, mi in mut2idx.items():
            df.append({
                'sample': s,
                'mutation': m,
                'fraction': spectra_fracs[si, mi],
                'haplotype': "B" if si in a_smp_i else "D" if si in b_smp_i else "H",
            })
    df = pd.DataFrame(df)
    sns.stripplot(
        df,
        x="mutation",
        y="fraction",
        hue="haplotype",
        ax=ax,
        palette="colorblind",
        dodge=True,
    )
    f.savefig('heatmap.png')
