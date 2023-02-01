import pandas as pd
import argparse
from scripts.run_ihd_scan import compute_spectra
import numpy as np
from compare_mutation_spectra import mutation_comparison
import matplotlib.pyplot as plt 
import seaborn as sns 
from sklearn.metrics.pairwise import cosine_distances, cosine_similarity
from statsmodels.formula.api import ols
import statsmodels.api as sm

PROJDIR = "/scratch/ucgd/lustre-work/quinlan/u1006375/proj-mutator-mapping"

geno = pd.read_csv(f"{PROJDIR}/data/genotypes/bxd.geno").set_index("marker")

markers = np.array(["rs30499894", "rs52263933"])
k = 1

geno_at_marker = geno.loc[markers].to_dict()

singletons = pd.read_csv(f"{PROJDIR}/data/singletons/bxd/annotated_filtered_singletons.csv")
singletons['count'] = 1

smp2generations = dict(zip(singletons['Strain'], singletons['n_generations']))

samples, mutations, spectra = compute_spectra(singletons, k=k)
smp2idx = dict(zip(samples, range(len(samples))))
mut2idx = dict(zip(mutations, range(len(mutations))))
samples_shared = list(set(samples).intersection(set(geno.columns[1:])))

smp2genotype = {}
for s in samples_shared:
    sorted_markers = sorted(geno_at_marker[s].items(), key=lambda t: t[0])
    #print (sorted_markers)
    marker_genos = "-".join([v for k,v in sorted_markers])
    smp2genotype[s] = marker_genos

a_smp_i = np.array([smp2idx[s] for s,v in smp2genotype.items() if v.split('-')[0] == "B"])
b_smp_i = np.array([smp2idx[s] for s,v in smp2genotype.items() if v.split('-')[0] == "D"])

a_spectra_sum = np.sum(spectra[a_smp_i], axis=0)
b_spectra_sum = np.sum(spectra[b_smp_i], axis=0)

if k == 3:
    mutation_comparison(b_spectra_sum, a_spectra_sum, mut2idx, vmin=-1, vmax=1)
elif k == 1:
    
    spectra_fracs = spectra / np.sum(spectra, axis=1)[:, np.newaxis]
    df = []
    for s in samples_shared:
        si = smp2idx[s]
        smp_geno = smp2genotype[s]
        smp_gens = smp2generations[s]
        if smp_gens < 0: continue
        if "H" in smp_geno: continue
        for m, mi in mut2idx.items():
            df.append({
                'sample': s,
                'Mutation type': m,
                'Fraction': spectra_fracs[si, mi],
                'Rate (per bp, per gen)': spectra[si, mi] / smp2generations[s] / 2.5e9,
                'Haplotype': smp_geno,
            })

    df = pd.DataFrame(df)
    print (df)

    lm = ols('Fraction ~ C(Haplotype)', data=df[df["Mutation type"] == "C>A"]).fit()
    table = sm.stats.anova_lm(lm, typ=2)
    ssq = table['sum_sq'].values
    print ((ssq[0] / np.sum(ssq)) * 100)

    f, (ax1, ax2) = plt.subplots(2, 1, figsize=(8, 12))
    for ax, kind in zip((ax1, ax2), ('Fraction', 'Rate (per bp, per gen)')):
        sns.boxplot(
            df,
            x="Mutation type",
            y=kind,
            hue="Haplotype",
            ax=ax,
            color="white",
            dodge=True,
            fliersize=0,
        )
        sns.stripplot(
            df,
            x="Mutation type",
            y=kind,
            hue="Haplotype",
            ax=ax,
            palette="colorblind",
            dodge=True,
        )
        # Get the handles and labels. For this example it'll be 2 tuples
        # of length 4 each.
    ax1.legend([],[], frameon=False)
    handles, labels = ax2.get_legend_handles_labels()
    l = plt.legend(handles[4:], labels[4:], title="Haplotypes at chr6 and chr4 peaks")
    f.tight_layout()
    f.savefig('heatmap.png', dpi=300)
