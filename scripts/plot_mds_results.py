import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import argparse
import numpy as np
from schema import IHDResultSchema, MarkerMetadataSchema
import scipy.stats as ss

plt.rc("font", size=16)

def main(args):

    # read in results of IHD scan and validate with pandera
    results = pd.read_csv(args.results)
    #IHDResultSchema.validate(results)

    # do the same with the genotype marker metadata
    markers = pd.read_csv(args.markers, dtype={"marker": str, "chromosome": str})
    MarkerMetadataSchema.validate(markers)

    results_merged = results.merge(markers, on="marker")

    chr_pref = any(["chr" in c for c in results_merged["chromosome"].unique()])
    if chr_pref:
        results_merged["chromosome"] = results_merged['chromosome'].apply(lambda c: int(c.split("chr")[1]))
    results_merged = results_merged[results_merged["chromosome"] != "X"]

    # get significant markers
    

    chrom_order = list(map(str, range(1, 23)))
    chrom_order.append("X")
    chrom_order_idx = dict(zip(chrom_order, range(len(chrom_order))))
    results_merged["sort_idx"] = results_merged["chromosome"].apply(lambda c: chrom_order_idx[c])
    results_merged.sort_values("sort_idx", inplace=True)

    #signif = results_merged.query("Distance >= 2.5 and genotype == 2")
    #signif.to_csv(f"{args.outpref}.significant_markers.csv", index=False)

    # plot manhattan

    for genotype in (0, 2):

        results_merged_sub = results_merged.query(f'genotype == {genotype}')
        results_merged_sub["MDS1_distance"] = ss.zscore(results_merged_sub["MDS1"].values)
        results_merged_sub["MDS2_distance"] = ss.zscore(results_merged_sub["MDS2"].values)

        f, (ax1, ax2) = plt.subplots(2, figsize=(14, 12))

        colors = ["cornflowerblue", "coral"]

        for d, ax in zip(("MDS1_distance", "MDS2_distance"), (ax1, ax2)):
            previous_max = 0
            xtick_positions, xticks = [], []

            for i, (
                    chrom,
                    chrom_df,
            ) in enumerate(results_merged_sub.groupby("chromosome", sort=False)):

                color_idx = i % 2
                xvals = chrom_df[args.colname].values + previous_max
                yvals = chrom_df[d].values
                ax.scatter(
                    xvals,
                    yvals,
                    s=50,
                    c=colors[color_idx],
                    ec="w",
                    lw=0.5,
                )

                previous_max += max(chrom_df[args.colname].values)
                xtick_positions.append(np.median(xvals))
                xticks.append(chrom)


        # for label, level, color, style in label_meta:
        #     max_dist = results_merged[f'{level}_percentile'].unique()[0]
        #     ax.axhline(y=max_dist, ls=style, c=color, label=label, lw=1.5)

                ax.set_xticks(xtick_positions)
                ax.set_xticklabels(xticks)
                ax.set_xlabel("Chromosome")
                ax.set_ylabel("Distance")
        #ax.legend()
        f.tight_layout()
        f.savefig(f"{args.outpref}.{genotype}.manhattan_plot.png", dpi=300)

    


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument(
        "--results",
        help=
        "Output of IHD scan, with an entry for every genotyped marker, in CSV format.",
    )
    p.add_argument(
        "--markers",
        help=
        "File containing metadata (cM position, Mb position, or both) about each genotyped marker.",
    )
    p.add_argument(
        "--outpref",
        type=str,
        default=".",
        help=
        "Prefix that will be prepended to output files and plots. Default is '.'",
    )
    p.add_argument(
        "-colname",
        default="Mb",
        help=
        "Name of the column in `--markers` that denotes the position you which to plot on the x-axis of the Manhattan plot. Default is 'Mb'",
    )
    p.add_argument(
        "-chrom",
        default=None,
        help="Chromosome to display separately in its own Manhattan plot.",
    )
    args = p.parse_args()
    main(args)