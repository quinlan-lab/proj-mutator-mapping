import pandas as pd
import argparse
import re
from collections import Counter
import json

def find_haplotype(genos: pd.DataFrame, sample: str) -> str:
    """
    figure out whether each strain has a B or D haplotype,
    or is heterozygous, at the genotype marker at the peak
    of the QTL on chromosome 4
    """

    genos_in_smp = genos[sample].values
    geno_freq = Counter(genos_in_smp)
    total = sum([i[1] for i in geno_freq.items()])
    most_freq_geno = "H"
    for g in ["B", "D"]:
        if geno_freq[g] > (total * 0.5): most_freq_geno = g[0]
        else: continue

    return most_freq_geno

def get_generation(gen: str) -> int:
    """
    given a Jackson Labs-formatted string that designates
    the inbreeding/backcrossing history of a strain, calculate
    the total number of generations of inbreeding a strain has
    undergone, and make a note of strains that have been backcrossed
    """

    # split generation designations by "+" symbols, which
    # indicate transitions between original facilities and JAX
    split = None
    try:
        split = re.split("(\d+)", gen)
    except TypeError:
        return 'NA'

    cur_gen = 0

    for i, e in enumerate(split):
        # add each number of filial generations to the
        # cumulative sum of generations
        if 'F' in e:
            cur_gen += int(split[i + 1])
        # "N"s in JAX designations indicate backcrossing
        # generations. we don't consider strains that have
        # been backcrossed, so return "NA"
        elif 'N' in e:
            return -1
        else:
            continue

    return int(cur_gen)

def main(args):
    combined = []
    for fh in args.singletons:
        df = pd.read_csv(fh)
        combined.append(df)

    combined = pd.concat(combined)

    # read in JSON file with file paths
    config_dict = None
    with open(args.config, "rb") as config:
        config_dict = json.load(config)

    geno_raw = pd.read_csv(config_dict['geno'])
    marker_info = pd.read_csv(config_dict['markers'])

    # merge genotype and marker information
    geno = geno_raw.merge(marker_info, on="marker")

    # get genotypes at top marker at chr4 eQTL
    rsids = ["rs52263933"]
    genos_at_markers = geno[geno['marker'].isin(rsids)]

    metadata = pd.read_excel(args.metadata)
    metadata['n_generations'] = metadata['Generation at sequencing'].apply(lambda g: get_generation(g))
    metadata = metadata[[
        "GeneNetwork name",
        "bam_name",
        "true_epoch",
        "n_generations",
    ]].dropna()
    metadata = metadata[metadata['n_generations'] != "NA"].astype({'n_generations': int})

    metadata = metadata.query('n_generations >= 20')
    #metadata = metadata[metadata['true_epoch'].isin([4])]

    combined_merged = combined.merge(metadata, left_on="bxd_strain", right_on="bam_name")
    combined_merged['Strain'] = combined_merged['GeneNetwork name']

    combined_merged['haplotype_at_qtl'] = combined_merged['Strain'].apply(
    lambda s: find_haplotype(genos_at_markers, s)
    if s in genos_at_markers.columns else "NA")

    #combined_merged = combined_merged[combined_merged['haplotype_at_qtl'] == "D"]

    combined_merged = combined_merged[combined_merged['Strain'] != "BXD68"]
    combined_merged.to_csv(args.out, index=False)

if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--singletons", nargs="+")
    p.add_argument("--metadata")
    p.add_argument("--config")
    p.add_argument("--out")
    args = p.parse_args()
    main(args)