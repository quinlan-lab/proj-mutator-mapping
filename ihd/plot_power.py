import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

plt.rc("font", size=12)

PROJDIR = "/Users/tomsasani/quinlanlab/proj-mutator-mapping"

df = pd.read_csv(f"{PROJDIR}/results.csv")
df['Power'] = df['pval'].apply(lambda p: p <= 0.05)
df["mutation_type"] = df["mutation_type"].apply(lambda m: m.replace(">", r"$\to$"))

replace_dict = {
    "mutation_type": "Mutation type",
    "n_haplotypes": "# haplotypes",
    "effect_size": "Mutator effect size",
    "n_mutations": "# mutations",
    "n_markers": "# markers",
    "tag_strength": "Tag strength",
}
df.rename(columns=replace_dict, inplace=True)


g = sns.FacetGrid(data=df, row="Tag strength", col="# haplotypes", aspect=1.5)
g.map(sns.lineplot, "Mutator effect size", "Power", "# mutations")
g.add_legend(title = "# of mutations\nper haplotype")
g.tight_layout()


g.savefig('sims.png', dpi=300)
