import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

plt.rc("font", size=16)

PROJDIR = "/Users/tomsasani/quinlanlab/proj-mutator-mapping"

markers = pd.read_csv(f"{PROJDIR}/data/genotypes/bxd.markers").query('chromosome != "X"')
max_length = markers.groupby('chromosome').agg({"Mb": max}).reset_index()
total_length = np.sum(max_length["Mb"])

max_length["percentage"] = max_length["Mb"] / total_length 
chrom2pct = dict(sorted(zip(max_length['chromosome'].astype(int), max_length['percentage']), key=lambda t: t[0], reverse=False))

df = pd.read_csv(f"{PROJDIR}/csv/heritability.csv").astype({"chromosome": int})

df["Mutation type"] = df["mutation"].apply(lambda m: m.replace('_', r"$\rightarrow$"))

palette = sns.color_palette("colorblind", 6)

ind = np.cumsum([v for k,v in chrom2pct.items()])
print (ind)
f, ax = plt.subplots(figsize=(12, 6))
for i, (mutation, mut_df) in enumerate(df.groupby("Mutation type")):
    mut_df = mut_df.sort_values("chromosome")
    cum_heritability = mut_df["heritability"].values[::-1]
    lw = 3 if mutation == "C" + r"$\rightarrow$" + "A" else 2 
    ls = "-" if mutation == "C" + r"$\rightarrow$" + "A" else "--"
    ax.plot(ind, cum_heritability, label=mutation, lw=lw, ls=ls, c=palette[i])

ax.legend(frameon=False)
ax.set_xticks(ind)
ax.set_xticklabels(k for k,v in chrom2pct.items())
ax.set_ylabel("Fraction of heritability explained")
ax.set_xlabel("Chromosome")
ax.set_title("Heritability of germline mutation fractions in the BXDs\n(after conditioning on the chr4 QTL)")
sns.despine(top=True, right=True, ax=ax)
f.tight_layout()
f.savefig('heritability.chr4_conditioned.png', dpi=300)
