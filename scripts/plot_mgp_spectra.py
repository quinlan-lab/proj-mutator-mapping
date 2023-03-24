import matplotlib.pyplot as plt 
import seaborn as sns 
import pandas as pd
import argparse

plt.rc('font', size=20)

def main(args):

    f, ax = plt.subplots(figsize=(14, 6))

    # read in dumont spectra dataframe
    dumont_tidy = pd.read_csv(args.spectra)

    # create color palette mapping
    palette = dict(
            zip(
                ["B-B", "B-D", "D-B", "D-D"],
                ["#398D84", "#E67F3A", "#EBBC2C", "#2F294A"],
            ))
    
    sns.boxplot(
        data=dumont_tidy.sort_values(["Mutation type", "Haplotypes"], ascending=True),
        x="Mutation type",
        y="Fraction",
        hue="Haplotypes",
        ax=ax,
        color="w",
        fliersize=0,
        linewidth=1,
    )
    sns.stripplot(
        data=dumont_tidy.sort_values(["Mutation type", "Haplotypes"], ascending=True),
        x="Mutation type",
        y="Fraction",
        hue="Haplotypes",
        ax=ax,
        palette=palette,
        dodge=True,
        ec="k",
        linewidth=1,
    )

    ax.set_ylabel('')#'Mutation rate (per bp, per gen.)')
    ax.set_ylim(0, 0.6)
    
    sns.despine(ax=ax, top=True, right=True)

    handles, labels = ax.get_legend_handles_labels()
    l = plt.legend(
        handles[4:],
        labels[4:],
        title="Genotypes at chr4 and chr6 mutator loci",
        frameon=False,
    )

    # change all spine widths
    for axis in ['top','bottom','left','right']:
        ax.spines[axis].set_linewidth(1.5)

    # increase tick width
    ax.tick_params(width=1.5)
    f.tight_layout()
    f.savefig(args.out, dpi=300)

if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument(
        "--spectra",
        help="""aggregate mutation spectra data in the MGP strains""",
    )
    p.add_argument(
        "--out",
        help="""name of output plot""",
    )
    args = p.parse_args()
    main(args)