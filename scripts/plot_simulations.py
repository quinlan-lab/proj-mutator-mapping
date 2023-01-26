import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

PROJDIR = "/scratch/ucgd/lustre-work/quinlan/u1006375/proj-mutator-mapping"

df = pd.read_csv(f"{PROJDIR}/results.csv")
df['pass'] = df['pval'].apply(lambda p: p <= 0.05)

group_cols = """n_haplotypes,frac_with_eQTL,augment_factor,mutation_count,mutation,n_markers""".split(",")
grouped_df = df.groupby(group_cols).agg({'pass': [np.mean, np.std]}).reset_index()

rename_cols = group_cols + ["mean_pass", "std_pass"]
grouped_df.columns = rename_cols

g = sns.FacetGrid(df, row="mutation", col="n_haplotypes")
g.map(sns.lineplot, "augment_factor", "pass", "mutation_count")
g.add_legend()
g.tight_layout()


g.savefig('sims.png', dpi=300)
