import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import statsmodels.stats.multitest as mt
import glob
import scipy.stats as ss

PROJDIR = "/scratch/ucgd/lustre-work/quinlan/u1006375/proj-bxd"

results = pd.read_csv(f"{PROJDIR}/o")
markers = pd.read_csv(f"{PROJDIR}/data/genotypes/BXD.markers")

results_merged = results.merge(markers, on="marker")

# try out method for 95% Bayes CI
# rescale cosine distances to be between 0 and 1
for chrom, chrom_df in results_merged.groupby('chromosome'):
    min_dist, max_dist = min(chrom_df['distance']), max(chrom_df['distance'])
    chrom_df['norm_distance'] = chrom_df['distance'].apply(lambda d: (d - min_dist) / (max_dist - min_dist))
    if chrom == "14":
        f, ax = plt.subplots()
        #cumsum = np.cumsum(chrom_df['norm_distance'])
        chrom_df['cumsum'] = chrom_df['distance'] / np.sum(chrom_df['distance'])
        sns.displot(chrom_df, x="Mb", y="distance")
        f.savefig('distplot.png')

# # compute approx. genetic length of each chromosome
# max_length_per_chrom = results_merged.groupby('chromosome').agg({'cM': max}).reset_index().rename(columns={'cM': 'genetic_length'})
# chrom_lens = max_length_per_chrom['genetic_length']
# max_length_per_chrom['frac'] = chrom_lens / np.sum(chrom_lens)

# # compute genome-adjusted p-values 
# results_merged = results_merged.merge(max_length_per_chrom, on="chromosome")

# results_merged['adj_pval'] = results_merged.apply(lambda row: 1 - ((1 - row['pval']) ** row['frac']), axis=1)
# #results_merged['-log10(p)'] = results_merged['adj_pval'].apply(lambda p: )

signif = results_merged[results_merged["distance"] >= results_merged["suggestive_percentile"]]
signif.to_csv(f"{PROJDIR}/significant_markers.csv", index=False)

g = sns.FacetGrid(results_merged, row="chromosome", sharex=False, aspect=2.5, sharey=True)
g.map(sns.scatterplot,
    "Mb",
    "distance",
    "k",
    palette="colorblind",
)
axes = g.axes[0]
for label, color in zip(('suggestive', 'significant'), ("dodgerblue", "firebrick")):
    max_dist = results_merged[f'{label}_percentile'].unique()[0]
    g.map(plt.axhline, y=max_dist, ls=":", c=color, label=f"Genome-wide {label} threshold")
g.add_legend()
g.tight_layout()
g.savefig('manhattan.png', dpi=300)