import pandas as pd
import argparse
import re

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

    metadata = pd.read_excel(args.metadata)
    metadata['n_generations'] = metadata['Generation at sequencing'].apply(lambda g: get_generation(g))
    metadata = metadata[[
        "GeneNetwork name",
        "bam_name",
        "true_epoch",
        "n_generations",
    ]].dropna()
    metadata = metadata[metadata['n_generations'] != "NA"].astype({'n_generations': int})

    #metadata = metadata.query('n_generations >= 20')
    #metadata = metadata[metadata['true_epoch'].isin([1,2,4])]
    #gn2bam = dict(zip(metadata['GeneNetwork name'], metadata['bam_name']))
    #bam2gn = {v:k for k,v in gn2bam.items()}

    combined_merged = combined.merge(metadata, left_on="bxd_strain", right_on="bam_name")

    #combined['Strain'] = combined['bxd_strain'].apply(lambda s: bam2gn[s] if s in bam2gn else "NA")
    #combined = combined[combined['Strain'] != "NA"]
    combined_merged['Strain'] = combined_merged['GeneNetwork name']
    combined_merged.to_csv(args.out, index=False)

if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--singletons", nargs="+")
    p.add_argument("--metadata")
    p.add_argument("--out")
    args = p.parse_args()
    main(args)