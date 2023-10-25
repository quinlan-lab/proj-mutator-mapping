import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import argparse

plt.rc("font", size=12)

def main(args):
    qtl_power = pd.read_csv(args.qtl_power)
    qtl_power["Method"] = "QTL"

    ihd_power = pd.read_csv(args.ihd_power)
    ihd_power['Power'] = ihd_power['pval'].apply(lambda p: p <= 0.05)
    ihd_power["Method"] = "AMSD"
    # use the cosine distance method in IHD to compare against QTL
    ihd_power = ihd_power[ihd_power["distance_method"] == "cosine"]

    correction = 7
    ihd_power["bonferroni_corr"] = correction

    combined_power = pd.concat([qtl_power, ihd_power]).reset_index()


    # ensure that we're comparing the same groups of simulations
    # between IHD and QTL mapping
    combined_power = combined_power[
        (combined_power["mutation_type"].isin(["C_A", "C_T"])) &
        (combined_power["n_haplotypes"] == 100) &
        (combined_power["bonferroni_corr"] == correction)]


    combined_power["mutation_type"] = combined_power["mutation_type"].apply(lambda m: m.replace("_", r"$\to$"))
    combined_power["effect_size"] = combined_power["effect_size"] / 100.
        
    replace_dict = {
        "mutation_type": "Mutation type",
        "n_haplotypes": "# haplotypes",
        "effect_size": "Mutator effect size",
        "n_mutations": "# mutations",
        "n_markers": "# markers",
        "exp_af": "Mutator allele freq.",
    }
    combined_power.rename(columns=replace_dict, inplace=True)

    row_name, col_name = "Mutation type", "# mutations"

    g = sns.FacetGrid(data=combined_power, row=row_name, col=col_name, aspect=1.5)
    g.map(sns.lineplot, "Mutator effect size", "Power", "Method", palette="colorblind")
    g.add_legend(title = "Method")
    g.tight_layout()

    g.savefig(args.out, dpi=300)

if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--qtl_power", help="""Results of power simulations using QTL mapping""")
    p.add_argument("--ihd_power", help="""Results of power simulations using IHD""")
    p.add_argument("--out", help="""Name of output plot""")
    args = p.parse_args()
    main(args)
