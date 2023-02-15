import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

plt.rc("font", size=16)

PROJDIR = "/Users/tomsasani/quinlanlab/proj-mutator-mapping"

df = pd.read_csv(f"{PROJDIR}/csv/heritability.csv")
df["Mutation type"] = df["mutation"].apply(lambda m: m.replace('_', r"$\rightarrow$"))

f, ax = plt.subplots(figsize=(12, 6))
sns.barplot(
    data=df,
    x="chromosome",
    y="heritability",
    hue="Mutation type",
    palette="colorblind",
    ec="k",
    lw=0.5,
    ax=ax,
)
handles, labels = ax.get_legend_handles_labels()
l = plt.legend(
    handles,
    labels,
    frameon=False,
)
ax.set_ylabel("Estimated heritability")
ax.set_xlabel("Chromosome")
ax.set_title("Heritability of 1-mer mutation fractions explained\nby markers on individual chromosomes")
sns.despine(top=True, right=True, ax=ax)
f.tight_layout()
f.savefig('heritability.chr4_samples_only.png', dpi=300)
