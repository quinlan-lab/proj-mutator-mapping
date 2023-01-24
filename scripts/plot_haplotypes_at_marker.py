import pandas as pd
import argparse
from compute_distance import compute_spectra
import numpy as np 
from compare_mutation_spectra import mutation_comparison
import matplotlib.pyplot as plt
import seaborn as sns

PROJDIR = "/scratch/ucgd/lustre-work/quinlan/u1006375/proj-bxd"

geno = pd.read_csv(f"{PROJDIR}/data/genotypes/BXD.geno")
markers = pd.read_csv(f"{PROJDIR}/data/genotypes/BXD.markers")
geno_merged = geno.merge(markers, on="marker")

singletons = pd.read_csv(f"{PROJDIR}/data/singletons/bxd/annotated_filtered_singletons.csv")
singletons['count'] = 1 

print (singletons.groupby('Strain').size().shape)

samples, mutations, spectra = compute_spectra(singletons, k=1)
#smp2idx = dict(zip(samples, range(len(samples))))
mut2idx = dict(zip(mutations, range(len(mutations))))
samples_shared = list(set(samples).intersection(set(geno.columns[1:])))

cols2use = ["marker", "chromosome", "Mb"]
cols2use.extend(samples_shared)
geno_merged = geno_merged[cols2use]

marker = "rs31187020"

geno_merged['is_focal_marker'] = geno_merged['marker'].apply(lambda m: 1 if m == marker else 0)
chrom = geno_merged.query('is_focal_marker == 1')['chromosome'].unique()[0]
pos = geno_merged.query('is_focal_marker == 1')['Mb'].unique()[0]

slop = 5
geno_to_plot = geno_merged[(geno_merged['chromosome'] == chrom) & (geno_merged['Mb'] > pos - slop) & (geno_merged['Mb'] < pos + slop)]
geno_to_plot.drop(columns="is_focal_marker", inplace=True)
geno_to_plot = geno_to_plot[geno_to_plot.columns[3:]].T.reset_index().rename(columns={'index': "Strain"})

# sort haplotypes by specified mutation fraction 
mutation = "C>A"
mut_idx = mut2idx[mutation]

spectra_fracs = spectra / np.sum(spectra, axis=1)[:, np.newaxis]
sorted_spectra_idxs = np.argsort(spectra_fracs[:, mut_idx])

smp2frac = dict(zip(samples, spectra_fracs[:, mut_idx]))
print (smp2frac)
geno_to_plot['frac'] = geno_to_plot['Strain'].apply(lambda s: smp2frac[s])
geno_to_plot_sorted = geno_to_plot.sort_values('frac')
print (geno_to_plot_sorted)
replace_dict = {'B': 0, 'D': 1, 'H': -1}
geno_to_plot_sorted.replace(to_replace=replace_dict, inplace=True)
geno_to_plot_sorted.fillna(-1, inplace=True)
geno_vals = geno_to_plot_sorted.values[:, 1:-1].astype(np.int8)
f, ax = plt.subplots()
sns.heatmap(geno_vals, ax=ax)
f.savefig('o.png')





