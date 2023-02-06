import pandas as pd
from run_ihd_scan import compute_spectra
import numpy as np
from compare_mutation_spectra import mutation_comparison
import matplotlib.pyplot as plt
import seaborn as sns
from statsmodels.formula.api import ols
import statsmodels.api as sm
from statannotations.Annotator import Annotator


plt.rc("font", size=14)

PROJDIR = "/Users/tomsasani/quinlanlab/proj-mutator-mapping"

geno = pd.read_csv(f"{PROJDIR}/data/genotypes/bxd.geno").set_index("marker")

markers = np.array(["rs46276051", "rs52263933"]) # chr6 and chr4
# markers = np.array(["rs6228198"]) # chr5
# markers = np.array(["rs13480950"]) # chr11
# markers = np.array(["rs30374203"])

k = 1

geno_at_marker = geno.loc[markers].to_dict()

mutations = pd.read_csv(f"{PROJDIR}/data/mutations/bxd/annotated_filtered_singletons.csv")

smp2generations = dict(zip(mutations['sample'], mutations['n_generations']))

samples, mutation_types, spectra = compute_spectra(mutations, k=k)
smp2idx = dict(zip(samples, range(len(samples))))
mut2idx = dict(zip(mutation_types, range(len(mutation_types))))
samples_shared = list(set(samples).intersection(set(geno.columns[1:])))

smp2genotype = {}
for s in samples_shared:
    sorted_markers = sorted(geno_at_marker[s].items(), key=lambda t: t[0])[::-1]
    marker_genos = "-".join([v for k,v in sorted_markers])
    smp2genotype[s] = marker_genos

a_smp_i = np.array([smp2idx[s] for s,v in smp2genotype.items() if v.split('-')[1] == "B"])
b_smp_i = np.array([smp2idx[s] for s,v in smp2genotype.items() if v.split('-')[1] == "D"])

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
                'Mutation type': m,#r"$\rightarrow$".join(m.split('>')),
                'Fraction': spectra_fracs[si, mi],
                'Rate': spectra[si, mi] / smp2generations[s] / 2.5e9,
                'Haplotype': smp_geno,
            })

    df = pd.DataFrame(df)

    # lm = ols("Rate ~ C(Haplotype, Treatment(reference='D-B'))", data=df[df["Mutation type"] == "C" + r"$\rightarrow$" + "A"]).fit()
    # print (lm.summary())
    # table = sm.stats.anova_lm(lm, typ=2)
    # ssq = table['sum_sq'].values
    # print ((ssq[0] / np.sum(ssq)) * 100)

    f, ax = plt.subplots(figsize=(10, 6))
    sns.boxplot(
        data=df,
        x="Mutation type",
        y="Rate",
        hue="Haplotype",
        ax=ax,
        color="white",
        dodge=True,
        fliersize=0,
    )
    sns.stripplot(
        data=df,
        x="Mutation type",
        y="Rate",
        hue="Haplotype",
        ax=ax,
        palette="colorblind",
        dodge=True,
    )

    pairs = [
        (("C>A", "B-B"), ("C>A", "B-D")),
        (("C>A", "B-D"), ("C>A", "D-D")),
        (("C>A", "D-B"), ("C>A", "D-D")),
    ]

    annotator = Annotator(
        ax,
        pairs,
        data=df,
        x="Mutation type",
        y="Rate",
        hue="Haplotype",
    )
    annotator.configure(
        #test='t-test_welch',
        test ="Mann-Whitney",
        comparisons_correction="BH",
        text_format='star',
    )
    annotator.apply_and_annotate()

    # Get the handles and labels. For this example it'll be 2 tuples
    # of length 4 each.
    handles, labels = ax.get_legend_handles_labels()
    l = plt.legend(
        handles[4:],
        labels[4:],
        title="Haplotypes at chr4 and chr6 peaks",
        frameon=False,
    )
    sns.despine(ax=ax, top=True, right=True)
    f.tight_layout()
    f.savefig('heatmap.png', dpi=300)
