import pandas as pd
from run_ihd_scan import compute_spectra
import numpy as np
from compare_mutation_spectra import mutation_comparison
import matplotlib.pyplot as plt
import seaborn as sns
from statsmodels.formula.api import ols
import statsmodels.api as sm
from statannotations.Annotator import Annotator
import scipy.stats as ss
from statsmodels.stats.multitest import multipletests


plt.rc("font", size=20)

PROJDIR = "/Users/tomsasani/quinlanlab/proj-mutator-mapping"

callable_kmers = pd.read_csv(f"{PROJDIR}/data/coverage/combined.callable_kmer_content.csv")

print (callable_kmers)

geno = pd.read_csv(f"{PROJDIR}/data/genotypes/bxd.geno").set_index("marker")

markers = np.array(["rs46276051", "rs52263933"]) # chr6 and chr4
#markers = np.array(["rs46276051"]) # chr6 and chr4

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


if k == 3:
    a_smp_i = np.array([smp2idx[s] for s,v in smp2genotype.items() if v.split('-')[1] == "B"])
    b_smp_i = np.array([smp2idx[s] for s,v in smp2genotype.items() if v.split('-')[1] == "D"])

    a_spectra_sum = np.sum(spectra[a_smp_i], axis=0)
    b_spectra_sum = np.sum(spectra[b_smp_i], axis=0)
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

        #rate_per_gen = spectra[si, mi] / smp2generations[s]
        callable_bp = callable_kmers[callable_kmers["GeneNetwork name"] == s]


        for m, mi in mut2idx.items():
            base_nuc = m.split(">")[0]
            n_callable_bp = callable_bp[callable_bp["nucleotide"] == base_nuc]["count"].values[0]
            rate = spectra[si, mi] / smp2generations[s] / n_callable_bp
            df.append({
                'sample': s,
                'Mutation type': r"$\rightarrow$".join(m.split('>')),
                'Count': spectra[si, mi],
                'Total': np.sum(spectra, axis=1)[si],
                'Fraction': spectra_fracs[si, mi],
                'Rate': rate,
                'Haplotype': smp_geno,
            })

    df = pd.DataFrame(df)
    df_grouped = df.drop(columns=["sample"]).groupby(["Haplotype", "Mutation type"]).agg(sum).reset_index()

    palette = ["#398D84", "#E67F3A", "#EBBC2C", "#2F294A"]

    # compare mutation rates using Chi-square tests
    comparisons = [
        ("B-B", "B-D"), # ogg1 against wt
        #("B-D", "D-D"), 
        ("D-B", "D-D"), # both against mutyh
        #("B-B", "D-B"),
    ]

    #comparisons = [("B", "D")]

    f, ax = plt.subplots(figsize=(10, 6))
    sns.boxplot(
        data=df,#[df["Haplotype"].isin(["D-B", "D-D"])],
        x="Mutation type",
        y="Rate",
        hue="Haplotype",
        ax=ax,
        color="white",
        fliersize=0,
    )
    sns.stripplot(
        data=df,#[df["Haplotype"].isin(["D-B", "D-D"])],
        x="Mutation type",
        y="Rate",
        palette=palette,
        hue="Haplotype",
        dodge=True,
        ax=ax,
    )
    sns.despine(ax=ax, top=True, right=True)
    for axis in ['top','bottom','left','right']:
        ax.spines[axis].set_linewidth(1.5)

    # increase tick width
    ax.tick_params(width=1.5)
    handles, labels = ax.get_legend_handles_labels()
    l = plt.legend(
        handles[4:],
        labels[4:],
        title="Genotypes at chr4 and chr6 peaks",
        frameon=False,
    )
    f.tight_layout()
    f.savefig("aggregate_mutation_spectra.jitter.png", dpi=300)


    pairs = []
    pvalues = []

    res = []
    for comparison_mut in df_grouped["Mutation type"].unique():

        for a_hap, b_hap in comparisons:
            a_hap_df = df_grouped[df_grouped["Haplotype"] == a_hap]
            b_hap_df = df_grouped[df_grouped["Haplotype"] == b_hap]
            # is focal mut and A hap
            a_fore = a_hap_df[a_hap_df["Mutation type"] == comparison_mut]["Count"].values[0]
            # is focal mut and B hap
            b_fore = b_hap_df[b_hap_df["Mutation type"] == comparison_mut]["Count"].values[0]
            # is not focal mut and A hap
            a_back = a_hap_df[a_hap_df["Mutation type"] != comparison_mut]["Count"].values.sum()
            # is not focal mut and B hap
            b_back = b_hap_df[b_hap_df["Mutation type"] != comparison_mut]["Count"].values.sum()

            _, p, _, _ = ss.chi2_contingency([
                [a_fore, b_fore],
                [a_back, b_back],
            ])

            #if comparison_mut != r"$\rightarrow$".join(["C", "A"]): 
            if p >= 0.05: continue
            pvalues.append(p)
            pairs.append(((comparison_mut, a_hap), (comparison_mut, b_hap)))

            res.append({
                "a_hap": a_hap,
                "b_hap": b_hap,
                "mut": comparison_mut,
                "p": p,
            })

    res_df = pd.DataFrame(res)
    _, adj_p, _, _ = multipletests(res_df["p"].values, method="fdr_bh")
    res_df["adj_p"] = adj_p

    df_grouped["Aggregate fraction"] = df_grouped["Count"] / df_grouped["Total"]

    plt.figure(figsize=(8, 6))
    ax = sns.barplot(
        data=df_grouped,
        x="Mutation type",
        y="Aggregate fraction",
        hue="Haplotype",
        #ax=ax,
        palette=palette,
        ec='k',
        lw=2,
    )
    annotator = Annotator(
        ax,
        pairs,
        data=df_grouped,
        x="Mutation type",
        y="Rate",
        hue="Haplotype",
    )

    annotator.configure(
        test=None,
        test_short_name=r"$\chi^{2}$",
        #text_format="full",
        #text_format="simple",
        text_format="star",
        #comparisons_correction="BH",
    ).set_pvalues(pvalues=pvalues).annotate()

    handles, labels = ax.get_legend_handles_labels()
    l = plt.legend(
        handles,
        labels,
        title="Haplotypes at chr4 and chr6 peaks",
        frameon=False,
    )
    #sns.set_style('ticks')
    sns.despine(ax=ax, top=True, right=True)
    # change all spines
    for axis in ['top','bottom','left','right']:
        ax.spines[axis].set_linewidth(1.5)

    # increase tick width
    ax.tick_params(width=1.5)
    plt.tight_layout()
    plt.savefig("aggregate_spectra_comparison.png", dpi=300)
    plt.savefig("aggregate_spectra_comparison.eps")
