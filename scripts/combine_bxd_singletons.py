import pandas as pd
import argparse
import re
from collections import Counter
import json


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

    # read in the singleton files and combine into
    # a single dataframe
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

    metadata = pd.read_excel(args.metadata)
    metadata['n_generations'] = metadata['Generation at sequencing'].apply(lambda g: get_generation(g))
    metadata = metadata[[
        "GeneNetwork name",
        "bam_name",
        "true_epoch",
        "n_generations",
    ]].dropna()
    metadata = metadata[metadata['n_generations'] != "NA"].astype({'n_generations': int})

    # merge the singleton data with relevant metadata
    combined_merged = combined.merge(
        metadata,
        left_on="bxd_strain",
        right_on="bam_name",
    )
    combined_merged['sample'] = combined_merged['GeneNetwork name']

    # if desired, only output data for samples that have D genotypes
    # at the chr4 QTL we identified previously
    if args.condition_on_mutator:
        geno_at_marker = geno[geno['marker'] == "rs52263933"]
        smp2geno = geno_at_marker.to_dict(orient="list")
        filtered_samples = [s for s,g in smp2geno.items() if g[0] == "D"]
        combined_merged = combined_merged[combined_merged["sample"].isin(filtered_samples)]

    # remove strains that haven't been inbred for a sufficient number of generations
    # and remove the BXD68 outlier
    combined_merged = combined_merged[(combined_merged['n_generations'] > 20)
                                      & (combined_merged["sample"] != "BXD68")]
    combined_merged['count'] = 1
    combined_merged.to_csv(args.out, index=False)


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument(
        "--singletons",
        nargs="+",
        help=
        """paths to each of the files containing per-chromosome singleton mutation info from Sasani et al. 2022""",
    )
    p.add_argument(
        "--metadata",
        help="""path to Excel file containing metadata about each BXD line""",
    )
    p.add_argument(
        "--config",
        help=
        """path to config JSON file containing paths to genotype and marker files for BXDs.""",
    )
    p.add_argument(
        "--out",
        help=
        """name of output file in which we'll store aggregate singleton mutation data""",
    )
    p.add_argument(
        "-condition_on_mutator",
        action="store_true",
        help=
        """if specified, only generate a combined singleton file using BXDs with D alleles at the chromosome 4 mutator locus discovered in Sasani et al. 2022""",
    )
    args = p.parse_args()
    main(args)
