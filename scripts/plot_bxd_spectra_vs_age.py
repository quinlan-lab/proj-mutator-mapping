import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import argparse
import statsmodels.api as sm
import statsmodels


plt.rc("font", size=18)

def main(args):

    spectra_df = pd.read_csv(args.spectra)

    palette = dict(
        zip(
            ["B-B", "B-D", "D-B", "D-D"],
            ["#398D84", "#E67F3A", "#EBBC2C", "#2F294A"],
        ))

    f, ax = plt.subplots(figsize=(8, 6))

    spectra_df = spectra_df[spectra_df["Mutation type"] == args.mutation_type.replace("_", ">")]
    # reformat mutation type
    spectra_df["Mutation type"] = spectra_df["Mutation type"].apply(lambda m: m.replace(">", r"$\to$"))

    for hap, hap_df in spectra_df.groupby("Haplotypes"):
        hap_df = hap_df.sort_values("Generations", ascending=True)
        X = hap_df["Generations"].values
        y = hap_df[args.phenotype].values
       
        sm.add_constant(X)
        link_func = statsmodels.genmod.families.links.identity()
        model = sm.GLM(y, X, family=sm.families.Poisson(link=link_func))
        results = model.fit()
        predictions = results.get_prediction()
        df_predictions = predictions.summary_frame(alpha=0.05)
        df_predictions["X"] = X

        muts_at_100 = results.get_prediction(100).summary_frame(alpha=0.05)
        mean_, ci_lo, ci_hi = (
            muts_at_100["mean"].values[0],
            muts_at_100["mean_ci_lower"].values[0],
            muts_at_100["mean_ci_upper"].values[0],
        )
        print (f"After 100 generations of inbreeding, {hap} \
                BXDs have accumulated {mean_} C>A mutations ({ci_lo} - {ci_hi})."                                                                                 )

        ax.plot(
            df_predictions["X"],
            df_predictions["mean"],
            label=hap,
            color=palette[hap],
            lw=2,
        )
        ax.fill_between(
            df_predictions["X"],
            y1=df_predictions["mean_ci_lower"],
            y2=df_predictions["mean_ci_upper"],
            color=palette[hap],
            alpha=0.25,
        )
        ax.scatter(
            X,
            y,
            ec="k",
            c=palette[hap],
            s=50,
        )

    sns.despine(ax=ax, top=True, right=True)
    for axis in ['top','bottom','left','right']:
        ax.spines[axis].set_linewidth(1.5)

    # increase tick width
    ylabel = r"C$\to$A " + f"mutation {args.phenotype.lower()}"
    ax.set_ylabel(ylabel)
    ax.set_xlabel("Number of generations of inbreeding")
    ax.tick_params(width=1.5)
    ax.legend(title="Genotypes at chr4 and chr6 peaks", frameon=False)
    f.tight_layout()
    f.savefig(args.out, dpi=300)

if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument(
        "--spectra",
        help="""Path to tidy dataframe of mutation spectra in BXDs""",
    )
    p.add_argument(
        "--out",
        help="""name of output file with tidy mutation spectra""",
    )
    p.add_argument(
        "-mutation_type",
        help="""Mutation type to plot. Default is C_A.""",
        type=str,
        default="C_A",
    )
    p.add_argument(
        "-phenotype",
        help="phenotype to use for the plot (options are Fraction, Count or Rate). Default is Count.",
        default="Count",
    )
    args = p.parse_args()
    main(args)