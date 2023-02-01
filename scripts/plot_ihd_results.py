import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import argparse
from schema import IHDResultSchema, MarkerMetadataSchema

def main(args):

    # read in results of IHD scan and validate with
    # pandera
    results = pd.read_csv(args.results)
    IHDResultSchema.validate(results)

    # do the same with the genotype marker metadata
    markers = pd.read_csv(args.markers)
    MarkerMetadataSchema.validate(markers)

    results_merged = results.merge(markers, on="marker")

    # get significant markers
    signif = results_merged[results_merged["distance"] >= results_merged["significant_percentile"]]
    signif.to_csv(f"{args.outpref}.significant_markers.csv", index=False)

    # plot manhattan
    g = sns.FacetGrid(results_merged, row="chromosome", sharex=False, aspect=2.5, sharey=True)
    g.map(sns.scatterplot,
        args.colname,
        "distance",
    )
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
    args = p.parse_args()
    main(args)