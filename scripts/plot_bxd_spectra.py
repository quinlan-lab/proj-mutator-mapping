import pandas as pd
import numpy as np
from compare_mutation_spectra import mutation_comparison
import matplotlib.pyplot as plt
import seaborn as sns
import argparse 


plt.rc("font", size=20)

def main(args):

    spectra_df = pd.read_csv(args.spectra)

    if args.k == 3:

        grouped_spectra_df = spectra_df.groupby(["Mutation type", "Haplotypes"]).agg({"Count": sum}).reset_index()

        grouped_spectra_df = grouped_spectra_df.pivot(index="Haplotypes", columns="Mutation type")
        mutation_types = [c[1].replace(r"$\rightarrow$", ">") for c in grouped_spectra_df.columns]
        mut2idx = dict(zip(mutation_types, range(len(mutation_types))))

        a_spectra_sum = grouped_spectra_df.loc["D-B"].values
        b_spectra_sum = grouped_spectra_df.loc["D-D"].values
        mutation_comparison(b_spectra_sum, a_spectra_sum, mut2idx, outname=args.out)

    elif args.k == 1:

        palette = dict(
        zip(
            ["B-B", "B-D", "D-B", "D-D"],
            ["#398D84", "#E67F3A", "#EBBC2C", "#2F294A"],
        ))

        f, ax = plt.subplots(figsize=(14, 6))

        # reformat mutation type
        spectra_df["Mutation type"] = spectra_df["Mutation type"].apply(lambda m: m.replace(">", r"$\to$"))

        # sort the spectra dataframe by mutation type and haplotype combination
        spectra_df.sort_values(["Mutation type", "Haplotypes"], ascending=True, inplace=True)

        spectra_df["ec"] = spectra_df["Mutation type"].apply(lambda m: "w" if m == "C>A" else "k")

        sns.boxplot(
            data=spectra_df,
            x="Mutation type",
            y=args.phenotype,
            hue="Haplotypes",
            ax=ax,
            color="white",
            fliersize=0,
        )
        sns.stripplot(
            data=spectra_df,
            x="Mutation type",
            y=args.phenotype,
            palette=palette,
            ec="k",
            linewidth=0.75,
            hue="Haplotypes",
            dodge=True,
            ax=ax,
        )
        sns.despine(ax=ax, top=True, right=True)
        for axis in ['top','bottom','left','right']:
            ax.spines[axis].set_linewidth(1.5)
        # bold the C>A mutation type
        for label in ax.get_xticklabels():
            if label.get_text() == r"C$\to$A": 
                label.set_fontweight('bold') 

        # increase tick width
        ax.tick_params(width=1.5)

        handles, labels = ax.get_legend_handles_labels()
        l = plt.legend(
            handles[4:],
            labels[4:],
            title="Genotypes at chr4 and chr6 peaks",
            frameon=False,
        )
        f.tight_layout()
        f.savefig(args.out, dpi=300)

if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument(
        "--spectra",
        help="""tidy dataframe of mutation spectra in BXDs""",
    )
    p.add_argument(
        "--out",
        help="""name of output file with tidy mutation spectra""",
    )
    p.add_argument(
        "-phenotype",
        help="phenotype to use for the plot (options are Fraction or Rate). Default is Fraction,",
        default="Fraction",
    )
    p.add_argument(
        "-k",
        help="""kmer context in which to compute mutation spectra""",
        type=int,
        default=1,
    )
    args = p.parse_args()
    main(args)