import pandas as pd 
import numpy as np

def add_fake_kmer(row):
    ref, alt = row['ref'], row['alt']
    return f"N{ref}N>N{alt}N"

PROJDIR = "/scratch/ucgd/lustre-work/quinlan/u1006375/proj-bxd"

genotypes = pd.read_csv(f"{PROJDIR}/data/singletons/cemee/S1_WS220_CeMEEv2_markerSet1.csv.gz")

ac = np.sum(genotypes.values[:, 4:], axis=1)
singleton_idxs = np.where(ac == 1)[0]
singletons = genotypes.iloc[singleton_idxs]
singletons['kmer'] = singletons.apply(lambda row: add_fake_kmer(row), axis=1)

cols2keep = [c for c in singletons.columns if c not in ("chrom", "pos", "ref", "alt")]

singletons_tidy = singletons[cols2keep].melt(id_vars="kmer", var_name="Strain", value_name="has_alt")

singletons = singletons_tidy.groupby(['Strain', 'kmer']).agg({'has_alt': sum}).reset_index().rename(columns={'has_alt': 'count'}).query('count > 0')

singletons.to_csv(f"{PROJDIR}/data/singletons/cemee/annotated_filtered_singletons.csv", index=False)

