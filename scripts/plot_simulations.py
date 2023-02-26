import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

plt.rc("font", size=12)

PROJDIR = "/Users/tomsasani/quinlanlab/proj-mutator-mapping"

df = pd.read_csv(f"{PROJDIR}/results.csv")
df['Power'] = df['pval'].apply(lambda p: p <= 0.05)

group_cols = [
    "# of haplotypes",
    "% with mutator",
    "Mutator effect size",
    "# of mutations",
    "Mutation type",
    "# genotyped markers",
]

g = sns.FacetGrid(df, row="Mutation type", col="# of haplotypes", aspect=1.5)
g.map(sns.lineplot, "Mutator effect size", "Power", "# of mutations")
g.add_legend(title = "# of mutations\nper haplotype")
g.tight_layout()


g.savefig('sims.png', dpi=300)
