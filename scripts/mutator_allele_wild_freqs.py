from cyvcf2 import VCF
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

plt.rc('font', size=14)

PROJDIR = "/Users/tomsasani/quinlanlab/proj-mutator-mapping"
vcf = VCF(f"{PROJDIR}/data/wild.vcf.gz", gts012=True)
smp2idx = dict(zip(vcf.samples, range(len(vcf.samples))))
idx2smp = {v:k for k,v in smp2idx.items()}

candidate_mutators = [
    "chr6:113328509-113328510",
    "chr6:108075852-108075853",
    "chr6:108076364-108076365",
]

candidate_names = [
    "Ogg1 (p.Thr95Ala)",
    "Setmar (p.Leu103Phe)",
    "Setmar (p.Ser273Arg)",
]

res = []

for region, name in zip(candidate_mutators, candidate_names):
    for v in vcf(region):
        gts = v.gt_types

        for idx, smp in idx2smp.items():

            res.append({
                'region': region,
                'Mutation': name,
                'sample': smp,
                'genotype': gts[idx],
            })

res_df = pd.DataFrame(res)
print (res_df)
res_df['Subspecies'] = res_df['sample'].apply(lambda s: s.split('_')[0])
res_df_grouped = res_df.groupby(['Mutation', 'Subspecies']).agg({
    'genotype':
    lambda g: sum(g) / (len(g) * 2)
}).reset_index().rename(columns={'genotype': "Allele frequency"}, )

palette = ["#398D84", "#E67F3A", "#EBBC2C", "#2F294A"]

f, ax = plt.subplots(figsize=(10, 6))
sns.barplot(
    data=res_df_grouped,
    x="Mutation",
    y="Allele frequency",
    hue="Subspecies",
    ax=ax,
    palette=palette,
    ec='k',
    lw=2,
)
sns.set_style("ticks")
sns.despine(ax=ax, top=True, right=True)
f.savefig("wild.af.png", dpi=300)
