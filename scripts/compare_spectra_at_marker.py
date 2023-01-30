import pandas as pd
import argparse
from compute_distance import compute_spectra
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

a_smp_i = np.array([smp2idx[s] for s,v in smp2genotype.items() if v.split('-')[1] == "B"])
b_smp_i = np.array([smp2idx[s] for s,v in smp2genotype.items() if v.split('-')[1] == "D"])

#a_spectra, b_spectra = spectra[a_smp_i], spectra[b_smp_i]
#all_spectra = np.concatenate((a_spectra, b_spectra))
#all_spectra_fracs = all_spectra / np.sum(all_spectra, axis=1)[:, np.newaxis]

a_spectra_sum = np.sum(spectra[a_smp_i], axis=0)
b_spectra_sum = np.sum(spectra[b_smp_i], axis=0)

if k == 3:
    mutation_comparison(b_spectra_sum, a_spectra_sum, mut2idx)
elif k == 1:
    f, ax = plt.subplots()
    spectra_fracs = spectra / np.sum(spectra, axis=1)[:, np.newaxis]
    df = []
    for s in samples_shared:
        si = smp2idx[s]
        smp_geno = smp2genotype[s]
        if "H" in smp_geno: continue
        for m, mi in mut2idx.items():
            df.append({
                'sample': s,
                'mutation': m,
                'fraction': spectra_fracs[si, mi],
                'haplotype': smp_geno,
            })

    df = pd.DataFrame(df)

    lm = ols('fraction ~ C(haplotype)', data=df.query('mutation == "C>A"')).fit()
    table = sm.stats.anova_lm(lm, typ=2)
    ssq = table['sum_sq'].values
    print ((ssq[0] / np.sum(ssq)) * 100)

    sns.boxplot(
        df,
        x="mutation",
        y="fraction",
        hue="haplotype",
        ax=ax,
        color="white",
        dodge=True,
    )
    sns.stripplot(
        df,
        x="mutation",
        y="fraction",
        hue="haplotype",
        ax=ax,
        palette="colorblind",
        dodge=True,
    )
    # Get the handles and labels. For this example it'll be 2 tuples
    # of length 4 each.
    handles, labels = ax.get_legend_handles_labels()
    

    # When creating the legend, only use the first two elements
    # to effectively remove the last two.
    l = plt.legend(handles[4:], labels[4:], title="Haplotypes at chr6 and chr4 peaks")#, bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    f.tight_layout()
    f.savefig('heatmap.png', dpi=300)
