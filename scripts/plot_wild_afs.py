from cyvcf2 import VCF
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import argparse

plt.rc('font', size=16)

def main(args):

    vcf = VCF(args.vcf, gts012=True)
    smp2idx = dict(zip(vcf.samples, range(len(vcf.samples))))
    idx2smp = {v:k for k, v in smp2idx.items()}

    candidate_mutators = [
        "chr6:108075852-108075853",
        "chr6:108076364-108076365",
        "chr6:113328509-113328510",
        "chr6:115849643-115849644"
    ]

    candidate_names = [
        "Setmar\n(p.Leu103Phe)",
        "Setmar\n(p.Ser273Arg)",
        "Ogg1\n(p.Thr95Ala)",
        "Mbd4\n(p.Asp129Asn)"
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

    res_df['Subspecies'] = res_df['sample'].apply(lambda s: s.split('_')[0])
    res_df_grouped = res_df.groupby(['Mutation', 'Subspecies']).agg({
        'genotype':
        lambda g: sum(g) / (len(g) * 2), 'sample': lambda g: len(g),
    }).reset_index().rename(columns={'genotype': "Allele frequency", 'sample': 'count'}, )
    res_df_grouped["Subspecies"] = res_df_grouped.apply(lambda row: f"{row['Subspecies']} ({row['count']})", axis=1)

    palette = ["#398D84", "#E67F3A", "#EBBC2C", "#2F294A"]

    f, ax = plt.subplots(figsize=(10, 6))
    sns.barplot(
        data=res_df_grouped,
        x="Mutation",
        y="Allele frequency",
        hue="Subspecies",
        ax=ax,
        palette=palette,
        ec='w',
        lw=2,
    )
    # change all spines
    for axis in ['top','bottom','left','right']:
        ax.spines[axis].set_linewidth(1.5)

    # increase tick width
    ax.tick_params(width=1.5)
    ax.legend(frameon=False, loc="upper center", ncol=2, alignment="center")
    ax.set_xlabel("")
    f.tight_layout()
    sns.despine(ax=ax, top=True, right=True)
    f.savefig(args.out, dpi=300)

if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--vcf", help="""Path to VCF containing wild mouse variants from Harr et al.""")
    p.add_argument("--out", help="""Name of output plot.""")
    args = p.parse_args()
    main(args)