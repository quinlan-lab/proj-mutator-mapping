import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.patches import Rectangle
import numpy as np

res = pd.read_csv("o.csv")
#X_idxs = res[res["chromosome"] == "X"].index.values

threshold = res["threshold"].values[0]

vals = res.drop(columns=[
    "marker",
    "chromosome",
    "threshold",
    "cM",
    "Mb",
]).values



vals += vals.T

f, ax = plt.subplots(figsize=(8, 6))
ax.imshow(vals, cmap="seismic")
# sns.heatmap(
#     vals,
#     ax=ax1,
# )
vals[vals <= threshold] = np.nan
print (np.where(~np.isnan(vals)))

x, y = 0, 0

mids = []
for i, (chrom, chrom_df) in enumerate(res.groupby("chromosome")):
    if chrom == "X": continue
    size = chrom_df.shape[0]
    
    rect = Rectangle(
        (x, y),
        width=size,
        height=size,
        edgecolor="k",
        lw=2,
        facecolor="none",
    )
    ax.add_patch(rect)
    mids.append(np.mean([x, (x + size)]))
    x += size
    y += size
print (res["chromosome"].unique())

ax.set_xticks(mids)
ax.set_xticklabels([c for c in res["chromosome"].unique() if c != "X"])
ax.set_yticks(mids)
ax.set_yticklabels([c for c in res["chromosome"].unique() if c != "X"])
ax.set_xlabel("Chromosome")
ax.set_ylabel("Chromosome")

ax.set_title("Raw epistasis signals between pairs of markers")
#ax2.set_title("Chromosome-wide significant epistasis signals (p = 0.05)")
f.tight_layout()
f.savefig('o.png', dpi=200)
