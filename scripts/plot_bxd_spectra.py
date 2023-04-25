import pandas as pd
from compare_mutation_spectra import mutation_comparison
import matplotlib.pyplot as plt
import seaborn as sns
import argparse
import statsmodels.api as sm
from statannotations.Annotator import Annotator

plt.rc("font", size=18)

def main(args):

    spectra_df = pd.read_csv(args.spectra)

    group_col = "Haplotypes"

    if args.k == 3:

        grouped_spectra_df = spectra_df.groupby(
            ["Mutation type", group_col]).agg({
                "Count": sum
            }).reset_index()

        grouped_spectra_df = grouped_spectra_df.pivot(index=group_col,
                                                      columns="Mutation type")
        mutation_types = [
            c[1].replace(r"$\rightarrow$", ">")
            for c in grouped_spectra_df.columns
        ]
        mut2idx = dict(zip(mutation_types, range(len(mutation_types))))

        a_spectra_sum = grouped_spectra_df.loc["D-B"].values
        b_spectra_sum = grouped_spectra_df.loc["D-D"].values
        mutation_comparison(b_spectra_sum,
                            a_spectra_sum,
                            mut2idx,
                            outname=args.out)

    elif args.k == 1:

        palette = dict(
            zip(
                ["B-B", "B-D", "D-B", "D-D"],
                ["#398D84", "#E67F3A", "#EBBC2C", "#2F294A"],
            ))

        f, ax = plt.subplots(figsize=(14, 6))
        lw = 0.5
        size = 3.5
        # if desired, subset to a single mutation type
        if args.mutation_type is not None:
            spectra_df = spectra_df[spectra_df["Mutation type"] ==
                                    args.mutation_type.replace("_", ">")]
            # fit linear model to get variance explained
            Y = spectra_df[args.phenotype]
            X = spectra_df[group_col]
            X = sm.add_constant(X)
            f, ax = plt.subplots(figsize=(8, 6))
            lw *= 1.5
            size *= 2

        # reformat mutation type
        spectra_df["Mutation type"] = spectra_df["Mutation type"].apply(
            lambda m: m.replace(">", r"$\to$"))

        # sort the spectra dataframe by mutation type and haplotype combination
        spectra_df.sort_values(["Mutation type", group_col],
                               ascending=True,
                               inplace=True)

        xval, hue, dodge = "Mutation type", group_col, True
        if args.mutation_type is not None:
            xval, hue, dodge = group_col, None, False

        comparisons = [("B-D", "B-B"), ("D-D", "D-B"), ("D-D", "B-D")]

        sns.boxplot(
            data=spectra_df,
            x=xval,
            y=args.phenotype,
            hue=hue,
            ax=ax,
            color="white",
            fliersize=0,
        )
        sns.stripplot(
            data=spectra_df,
            x=xval,
            y=args.phenotype,
            palette=palette,
            edgecolor="black",
            linewidth=lw,
            hue=hue,
            size=size,
            dodge=dodge,
            ax=ax,
        )

        if args.mutation_type is not None:
            # annotate plot
            annotator = Annotator(
                    ax,
                    comparisons,
                    data=spectra_df,
                    x=xval,
                    y=args.phenotype,
                    hue=hue,
                )
            annotator.configure(
                test="Mann-Whitney",
                show_test_name=False,
                text_format="simple",
                text_offset=10,
                comparisons_correction=None,
            ).apply_and_annotate()

        sns.despine(ax=ax, top=True, right=True)
        for axis in ['top', 'bottom', 'left', 'right']:
            ax.spines[axis].set_linewidth(1.5)

        # increase tick width
        ax.tick_params(width=1.5)
        if args.mutation_type is not None:
            ax.set_ylabel(r"C$\to$A" + f" mutation {args.phenotype.lower()}")
            ax.set_xlabel("Genotypes at chr4 and chr6 peaks")

        handles, labels = ax.get_legend_handles_labels()
        handles, labels = (handles[4:], labels[4:])
        title = "Genotypes at chr4 and chr6 peaks"
        if args.mutation_type is not None:
            handles, labels, title = [], [], None
        l = plt.legend(
            handles, labels,
            title=title,
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
        help=
        "phenotype to use for the plot (options are Fraction or Rate). Default is Fraction,",
        default="Fraction",
    )
    p.add_argument(
        "-k",
        help="""kmer context in which to compute mutation spectra""",
        type=int,
        default=1,
    )
    p.add_argument(
        "-mutation_type",
        help=
        """If specified, only plot mutation fractions for a single mutation type.""",
        default=None,
    )
    args = p.parse_args()
    main(args)