import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import seaborn as sns
import argparse
import numpy as np
from schema import IHDResultSchema, MarkerMetadataSchema

plt.rc("font", size=18)

def main(args):

    # read in results of IHD scan and validate with pandera
    results = pd.read_csv(args.results)
    IHDResultSchema.validate(results)

    # do the same with the genotype marker metadata
    markers = pd.read_csv(args.markers, dtype={"marker": str, "chromosome": str})
    MarkerMetadataSchema.validate(markers)

    results_merged = results.merge(markers, on="marker").astype({"chromosome": str})

    chr_pref = any(["chr" in c for c in results_merged["chromosome"].unique()])
    if chr_pref:
        results_merged["chromosome"] = results_merged['chromosome'].apply(lambda c: int(c.split("chr")[1]))
    results_merged = results_merged[results_merged["chromosome"] != "X"]

    pctile_label = f"95th_percentile"

    # get significant markers
    signif = results_merged[results_merged["Distance"] >= results_merged[pctile_label]]
    signif.to_csv(f"{args.outpref}.significant_markers.csv", index=False)

    chrom_order = list(map(str, range(1, 23)))
    chrom_order.append("X")
    chrom_order_idx = dict(zip(chrom_order, range(len(chrom_order))))
    results_merged["sort_idx"] = results_merged["chromosome"].apply(lambda c: chrom_order_idx[c])
    results_merged.sort_values("sort_idx", inplace=True)

    results_merged["is_significant"] = results_merged["Distance"] >= results_merged[pctile_label]
    results_merged["ec"] = results_merged["is_significant"].apply(lambda s: "k" if s else "w")
    results_merged["lw"] = results_merged["is_significant"].apply(lambda s: 1 if s else 0.5)

    # plot manhattan
    f, ax = plt.subplots(figsize=(16, 6))

    max_threshold_dist = np.max(results_merged[pctile_label])

    ax.axhline(
        y=max_threshold_dist,
        ls=":",
        c="grey",
        label="Genome-wide significance threshold " +
        r"$\left(p \leq 0.05\right)$",
        lw=2,
    )

    previous_max = 0
    xtick_positions, xticks = [], []

    colors = ["#E6803C", "#398D84"]

    for i, (
            chrom,
            chrom_df,
    ) in enumerate(results_merged.groupby("chromosome", sort=False)):

        color_idx = i % 2
        xvals = chrom_df[args.colname].values + previous_max
        yvals = chrom_df["Distance"].values

        ax.scatter(
            xvals,
            yvals,
            s=75,
            c=colors[color_idx],
            ec=chrom_df["ec"].values,
            lw=chrom_df["lw"].values,
        )

        previous_max += max(chrom_df[args.colname].values)
        xtick_positions.append(np.median(xvals))
        xticks.append(chrom)


    ax.set_xticks(xtick_positions)
    ax.set_xticklabels(xticks)

    # change all spines
    for axis in ['top','bottom','left','right']:
        ax.spines[axis].set_linewidth(2.)

    # increase tick width
    ax.tick_params(width=2.)

    max_yval = max([
        max(results_merged["Distance"]) * 1.05,
        max_threshold_dist,
    ])

    min_yval = min([
        min(results_merged["Distance"]) * 1.05,
        max_threshold_dist,
    ])

    yticks = np.linspace(min_yval, max_yval, num=8)
    ytick_labels = [f"{x:.1e}" for x in yticks]
    ytick_labels = [round(x, 1) for x in yticks]

    #ytick_labels[0] = "0"
    ax.set_yticks(yticks[1:])
    ax.set_yticklabels(ytick_labels[1:])
    #sns.set_style('ticks')
    sns.despine(ax=ax, top=True, right=True)
    ax.set_xlabel("Chromosome")
    ax.set_ylabel("Cosine distance residual")
    ax.legend(frameon=False)
    f.tight_layout()
    f.savefig(f"{args.outpref}.manhattan_plot.png", dpi=300)
    f.savefig(f"{args.outpref}.manhattan_plot.eps")

    if args.chrom is not None:
        f, ax = plt.subplots(figsize=(10, 5))

        results_merged_chr = results_merged[results_merged["chromosome"] == args.chrom]
        max_threshold_dist = np.max(results_merged_chr[pctile_label])

        ax.axhline(
            y=max_threshold_dist,
            ls=":",
            c="grey",
            label="Genome-wide significance threshold " +
            r"$\left(p \leq 0.05\right)$",
            lw=2,
        )

        xvals = results_merged_chr[args.colname].values
        yvals = results_merged_chr["Distance"].values

        ax.scatter(
            xvals,
            yvals,
            s=75,
            c=colors[0],
            ec=results_merged_chr["ec"].values,
            lw=results_merged_chr["lw"].values,
        )

        ax.legend()
        f.savefig(f"{args.outpref}.{args.chrom}.manhattan_plot.png", dpi=300)
        #f.savefig(f"{args.outpref}.{args.chrom}.manhattan_plot.eps")


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