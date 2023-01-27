import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import statsmodels.stats.multitest as mt
import glob
import scipy.stats as ss
import argparse
from sklearn.linear_model import LinearRegression

def main(args):

    results = pd.read_csv(args.results)
    markers = pd.read_csv(args.markers)

    results_merged = results.merge(markers, on="marker")

    # fit PCA on genotype matrix
    # X, y = results_merged["genetic_difference"].values, results_merged["distance"].values
    # reg = LinearRegression().fit(X.reshape(-1, 1), y)
    # preds = reg.predict(X.reshape(-1, 1))
    # resids = preds - y
    # results_merged["residuals"] = resids * -1

    # f, ax = plt.subplots()
    # ax.scatter(y, preds)
    # f.savefig('resids.png', dpi=300)


    # # compute approx. genetic length of each chromosome
    # max_length_per_chrom = results_merged.groupby('chromosome').agg({'cM': max}).reset_index().rename(columns={'cM': 'genetic_length'})
    # chrom_lens = max_length_per_chrom['genetic_length']
    # max_length_per_chrom['frac'] = chrom_lens / np.sum(chrom_lens)

    # # compute genome-adjusted p-values
    # results_merged = results_merged.merge(max_length_per_chrom, on="chromosome")

    # results_merged['adj_pval'] = results_merged.apply(lambda row: 1 - ((1 - row['pval']) ** row['frac']), axis=1)
    # results_merged['adjusted_sig'] = results_merged['adj_pval'].apply(lambda p: 1 if p <= 0.05 else 0)

    results_merged['odds_ratio'] = results_merged['distance'] / results_merged['null_mean']
    results_merged['suggestive_odds_ratio'] = results_merged['suggestive_percentile'] / results_merged['null_mean']
    results_merged['significant_odds_ratio'] = results_merged['significant_percentile'] / results_merged['null_mean']

    signif = results_merged[results_merged["distance"] >= results_merged["suggestive_percentile"]]
    signif.to_csv(f"{args.outpref}.significant_markers.csv", index=False)

    g = sns.FacetGrid(results_merged, row="chromosome", sharex=False, aspect=2.5, sharey=True)
    g.map(sns.scatterplot,
        "Mb",
        "distance",
        "k",
        palette="colorblind",
    )
    axes = g.axes[0]
    for label, level, color in zip(
        ('Suggestive distance threshold', 'Significant distance threshold'),
        ('suggestive', 'significant'),
        ("dodgerblue", "firebrick"),
    ):
        max_dist = results_merged[f'{level}_percentile'].unique()[0]
        g.map(plt.axhline, y=max_dist, ls=":", c=color, label=label)
    g.add_legend()
    g.tight_layout()
    g.savefig(f"{args.outpref}.manhattan_plot.png", dpi=300)

if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--results")
    p.add_argument("--markers")
    p.add_argument("--outpref", type=str, default=".")
    args = p.parse_args()
    main(args)