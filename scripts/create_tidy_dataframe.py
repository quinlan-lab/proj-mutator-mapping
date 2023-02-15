import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import argparse
import json
import sys
from utils import compute_spectra
from skbio.stats.composition import clr, ilr
from combine_bxd_singletons import find_haplotype

from schema import IHDResultSchema, MutationSchema


PROJDIR = "/Users/tomsasani/quinlanlab/proj-mutator-mapping"

def main():

    # read in JSON file with file paths
    config_dict = None
    with open(f"{PROJDIR}/data/json/bxd.json", "rb") as config:
        config_dict = json.load(config)

    # read in genotype info
    geno = pd.read_csv(config_dict['geno'])

    # read in singleton data and validate with pandera
    mutations = pd.read_csv(f"{PROJDIR}/data/mutations/bxd/annotated_filtered_singletons.csv")
    MutationSchema.validate(mutations)

    smp2generation = dict(zip(mutations["sample"], mutations["n_generations"]))

    callable_kmers = pd.read_csv(f"{PROJDIR}/data/coverage/combined.callable_kmer_content.csv")

    samples = mutations['sample'].unique()

    # for the null permutation test, shuffle the rows of the spectra
    # dataframe every time. otherwise keep it the same.
    samples, mutation_types, spectra = compute_spectra(mutations, k=1)

    spectra_sums = spectra

    spectra_df = []
    for si, s in enumerate(samples):
        s_meta = callable_kmers[callable_kmers["GeneNetwork name"] == s]
        n_generations = smp2generation[s]
        for mi, m in enumerate(mutation_types):
            sm_count = spectra[si, mi]
            base_nuc = m.split(">")[0]
            n_callable_bp = s_meta[s_meta["nucleotide"] == base_nuc]["count"].values[0]
            spectra_df.append({
                "sample": s,
                "mutation": m.replace(">", "_"),
                "count": sm_count,
                "rate": sm_count / n_generations / n_callable_bp,
            })
    spectra_df = pd.DataFrame(spectra_df)

    spectra_grouped = spectra_df.groupby(["sample", "mutation"]).agg(sum).reset_index()
    spectra_totals = spectra_grouped.groupby("sample").agg({"rate": sum}).reset_index().rename(columns={"rate": "total_rate"})
    spectra_merged = spectra_grouped.merge(spectra_totals, on="sample")#["Aggregate fraction"] = spectra_grouped["Count"] / spectra_grouped["Total"]

    spectra_merged["fraction"] = spectra_merged["rate"] / spectra_merged["total_rate"]

    spectra_wide = spectra_merged.pivot(
        index="sample",
        columns="mutation",
        values="fraction",
    ).reset_index()

    # add covariates
    for rsid, qtl_name in zip(["rs52263933", "rs31412077"], ["chr4", "chr6"]):
        # get genotypes at top marker at QTL
        genos_at_markers = geno[geno['marker'] == rsid]

        spectra_wide[f'haplotype_at_{qtl_name}_qtl'] = spectra_wide['sample'].apply(
            lambda s: find_haplotype(genos_at_markers, s)
            if s in genos_at_markers.columns else -1)

    geno_samples = geno.columns[1:]
    sample_overlap = list(set(geno_samples).intersection(set(spectra_wide["sample"].to_list())))
    spectra_wide = spectra_wide[spectra_wide['sample'].isin(sample_overlap)]

    spectra_wide = spectra_wide[spectra_wide["haplotype_at_chr4_qtl"] == 1]

    raw_spectra = spectra_wide[[c.replace(">", "_") for c in mutation_types]].values
    #clr_spectra = clr(raw_spectra)
    # spectra_wide_clr = pd.DataFrame(
    #     clr_spectra,
    #     columns=[c.replace(">", "_") for c in mutation_types],
    # )
    # spectra_wide_clr['sample'] = spectra_wide["sample"]
    # add covariates
    # for rsid, qtl_name in zip(["rs52263933", "rs31412077"], ["chr4", "chr6"]):
    #     # get genotypes at top marker at QTL
    #     genos_at_markers = geno[geno['marker'] == rsid]

    #     spectra_wide_clr[f'haplotype_at_{qtl_name}_qtl'] = spectra_wide_clr['sample'].apply(
    #         lambda s: find_haplotype(genos_at_markers, s)
    #         if s in genos_at_markers.columns else -1)


    spectra_wide.dropna().to_csv(f"{PROJDIR}/csv/tidy_spectra.csv", index=False)

    # prepare cross info for heritability as well
    cross_info = spectra_wide[["sample"]].rename(columns={"sample": "id"})
    cross_info["cross_direction"] = "BxD"
    cross_info.to_csv(f"{PROJDIR}/data/Rqtl_data/bxd.cross_info.csv", index=False)

if __name__ == "__main__":
    main()
