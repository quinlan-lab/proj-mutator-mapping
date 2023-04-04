import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import argparse

plt.rc("font", size=12)

def main(args):
    ihd_power = pd.read_csv(args.ihd_power)
    ihd_power['Power'] = ihd_power['pval'].apply(lambda p: p <= 0.05)

    ihd_power = ihd_power[(ihd_power["mutation_type"].isin(["TCC_TTC", "C_T"])) & (ihd_power["n_haplotypes"] == 100)]

    ihd_power["mutation_type"] = ihd_power["mutation_type"].apply(lambda m: m.replace("_", r"$\to$"))
    ihd_power["effect_size"] = ihd_power["effect_size"] / 100.

    replace_dict = {
        "mutation_type": "Mutation type",
        "n_haplotypes": "# haplotypes",
        "effect_size": "Mutator effect size",
        "n_mutations": "# mutations",
        "n_markers": "# markers",
        "tag_strength": "Tag strength",
        "distance_method": "Distance method",
    }
    ihd_power.rename(columns=replace_dict, inplace=True)

    g = sns.FacetGrid(data=ihd_power, row="# mutations", col="Mutation type", aspect=1.5)
    g.map(sns.lineplot, "Mutator effect size", "Power", "Distance method", palette="colorblind")
    g.add_legend(title = "Distance method")
    g.tight_layout()

    g.savefig(args.out, dpi=300)

if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--qtl_power", help="""Results of power simulations using QTL mapping""")
    p.add_argument("--ihd_power", help="""Results of power simulations using IHD""")
    p.add_argument("--out", help="""Name of output plot""")
    args = p.parse_args()
    main(args)
