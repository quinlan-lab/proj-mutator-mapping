import pandas as pd
import argparse
from scripts.run_ihd_scan import compute_spectra
import numpy as np 
from compare_mutation_spectra import mutation_comparison
import matplotlib.pyplot as plt
import seaborn as sns

PROJDIR = "/scratch/ucgd/lustre-work/quinlan/u1006375/proj-mutator-mapping"

geno = pd.read_csv(f"{PROJDIR}/data/genotypes/BXD.geno")
markers = pd.read_csv(f"{PROJDIR}/data/genotypes/BXD.markers")
geno_merged = geno.merge(markers, on="marker")

singletons = pd.read_csv(f"{PROJDIR}/data/singletons/bxd/annotated_filtered_singletons.csv")
singletons['count'] = 1

epoch2smp = singletons.drop_duplicates('Strain').groupby('true_epoch')['Strain'].apply(list).to_dict()

samples, mutations, spectra = compute_spectra(singletons, k=1)
mut2idx = dict(zip(mutations, range(len(mutations))))
samples_shared = list(set(samples).intersection(set(geno.columns[1:])))

cols2use = ["marker", "chromosome", "Mb"]
cols2use.extend(samples_shared)
geno_merged = geno_merged[cols2use]

geno_merged = geno_merged[geno_merged['chromosome'] == "4"]

replace_dict = {'B': 0., 'D': 1., 'H': np.nan}

geno_merged.replace(to_replace=replace_dict, inplace=True)

f, ax = plt.subplots()

for epoch, smps in epoch2smp.items():
    if len(smps) < 10: continue
    cols2use = ["marker", "chromosome", "Mb"]
    cols2use.extend([s for s in smps if s in samples_shared])
    geno_sub = geno_merged[cols2use]
    genotypes = geno_sub.values[:, 3:].astype(np.float32)
    ac = np.nansum(genotypes, axis=1)
    an = np.sum(~np.isnan(genotypes), axis=1)
    afs = ac / an
    positions = geno_sub["Mb"].values 
    ax.plot(positions, afs, label=epoch)
ax.legend()
f.tight_layout()
f.savefig('afs.png', dpi=300)








