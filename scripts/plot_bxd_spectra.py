import pandas as pd
import numpy as np
from compare_mutation_spectra import mutation_comparison
import matplotlib.pyplot as plt
import seaborn as sns


plt.rc("font", size=20)

def main(args):

    spectra_df = pd.read_csv(args.spectra)
    print (spectra_df)

    # if args.k == 3:
    #     a_smp_i = np.array([smp2idx[s] for s,v in smp2genotype.items() if v == "D-B"])
    #     b_smp_i = np.array([smp2idx[s] for s,v in smp2genotype.items() if v == "D-D"])


    #     a_spectra_sum = np.sum(spectra[a_smp_i], axis=0)
    #     b_spectra_sum = np.sum(spectra[b_smp_i], axis=0)
    #     mutation_comparison(b_spectra_sum, a_spectra_sum, mut2idx, outname="heatmap.png")#, vmin=-1, vmax=1, )
    if args.k == 1:

        palette = dict(
        zip(
            ["B-B", "B-D", "D-B", "D-D"],
            ["#398D84", "#E67F3A", "#EBBC2C", "#2F294A"],
        ))

        f, ax = plt.subplots(figsize=(14, 6))

        sns.boxplot(
            data=spectra_df.sort_values(["Mutation type", "Haplotypes"], ascending=True),
            x="Mutation type",
            y=args.phenotype,
            hue="Haplotypes",
            ax=ax,
            color="white",
            fliersize=0,
        )
        sns.stripplot(
            data=spectra_df.sort_values(["Mutation type", "Haplotypes"], ascending=True),
            x="Mutation type",
            y=args.phenotype,
            palette=palette,
            ec='k',
            linewidth=0.75,
            hue="Haplotypes",
            dodge=True,
            ax=ax,
        )
        sns.despine(ax=ax, top=True, right=True)
        for axis in ['top','bottom','left','right']:
            ax.spines[axis].set_linewidth(1.5)

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
