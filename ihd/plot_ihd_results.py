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
        results_merged["chromosome"] = results_merged["chromosome"].apply(
            lambda c: int(c.split("chr")[1])
        )

    results_merged = results_merged[~results_merged["chromosome"].isin(["Y", "M"])]

    pctile_label = f"95th_percentile"

    chrom_order = list(map(str, range(1, 23)))
    chrom_order.append("X")
    chrom_order_idx = dict(zip(chrom_order, range(len(chrom_order))))
    results_merged["sort_idx"] = results_merged["chromosome"].apply(
        lambda c: chrom_order_idx[c]
    )
    results_merged.sort_values(["sort_idx", args.colname], inplace=True)

    results_merged["is_significant"] = (
        results_merged["Distance"] >= results_merged[pctile_label]
    )
    results_merged["ec"] = results_merged["is_significant"].apply(
        lambda s: "k" if s else "w"
    )
    results_merged["lw"] = results_merged["is_significant"].apply(
        lambda s: 1 if s else 0.5
    )

    # plot manhattan
    f, ax = plt.subplots(figsize=(16, 6))

    max_threshold_dist = np.max(results_merged[pctile_label])

    ax.axhline(
        y=max_threshold_dist * args.scale,
        ls=":",
        c="grey",
        label="Genome-wide significance threshold " + r"$\left(p = 0.05\right)$",
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
        yvals = chrom_df["Distance"].values * args.scale

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

    n_yticks = 6
    max_ypos = np.max(
        [
            np.max(results_merged["Distance"].values * args.scale),
            max_threshold_dist,
        ]
    )
    min_ypos = np.min(
        [
            np.min(results_merged["Distance"].values * args.scale),
            max_threshold_dist,
        ]
    )

    ytick_pos = np.linspace(min_ypos * 0.95, max_ypos * 1.05, n_yticks)
    ytick_labs = [f"{yval:.1e}" for yval in ytick_pos]

    ax.set_yticks(ytick_pos)
    ax.set_yticklabels(ytick_labs)

    # change all spines
    for axis in ["top", "bottom", "left", "right"]:
        ax.spines[axis].set_linewidth(1.5)

    # increase tick width
    ax.tick_params(width=2.0)

    sns.despine(ax=ax, top=True, right=True)
    ax.set_xlabel("Chromosome")
    ax.set_ylabel(f"Adjusted cosine distance")
    ax.legend(frameon=False)
    f.tight_layout()
    f.savefig(args.out, dpi=300)

    if args.chrom is not None:
        f, ax = plt.subplots(figsize=(8, 5))

        results_merged_chr = results_merged[results_merged["chromosome"] == args.chrom]
        max_threshold_dist = np.max(results_merged_chr[pctile_label])

        ax.axhline(
            y=max_threshold_dist * args.scale,
            ls=":",
            c="grey",
            label="Genome-wide significance threshold " + r"$\left(p = 0.05\right)$",
            lw=2,
        )

        xvals = results_merged_chr[args.colname].values
        yvals = results_merged_chr["Distance"].values * args.scale

        ax.scatter(
            xvals,
            yvals,
            s=75,
            c=colors[1],
            ec=results_merged_chr["ec"].values,
            lw=results_merged_chr["lw"].values,
        )

        ax.set_xlabel(f"Position on chr{args.chrom} (Mbp)")
        ax.set_ylabel("Adjusted cosine distance")
        f.tight_layout()

        f.savefig(args.out, dpi=300)


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument(
        "--results",
        help="Output of IHD scan, with an entry for every genotyped marker, in CSV format.",
    )
    p.add_argument(
        "--markers",
        help="File containing metadata (cM position, Mb position, or both) about each genotyped marker.",
    )
    p.add_argument(
        "--out",
        type=str,
        default="ihd.png",
        help="Name of output plot. Default is 'ihd.png'",
    )
    p.add_argument(
        "-colname",
        default="Mb",
        help="Name of the column in `--markers` that denotes the position you which to plot on the x-axis of the Manhattan plot. Default is 'Mb'",
    )
    p.add_argument(
        "-chrom",
        default=None,
        help="Chromosome to display separately in its own Manhattan plot.",
    )
    p.add_argument(
        "-scale",
        type=int,
        default=1,
        help="""Scale the cosine distance values by the specified amount to make visualization a bit easier on the y-axis. Default is 1000.""",
    )
    args = p.parse_args()
    main(args)
