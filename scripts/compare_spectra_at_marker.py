import pandas as pd
import argparse
from compute_distance import compute_spectra
import numpy as np 
from compare_mutation_spectra import mutation_comparison

PROJDIR = "/scratch/ucgd/lustre-work/quinlan/u1006375/proj-bxd"

geno = pd.read_csv(f"{PROJDIR}/data/genotypes/BXD.geno").set_index("marker")

marker = "rs46045549"

geno_at_marker = geno.loc[marker]

singletons = pd.read_csv(f"{PROJDIR}/data/singletons/bxd/annotated_filtered_singletons.csv")
singletons['count'] = 1 

samples, mutations, spectra = compute_spectra(singletons, k=3)
smp2idx = dict(zip(samples, range(len(samples))))
mut2idx = dict(zip(mutations, range(len(mutations))))
samples_shared = list(set(samples).intersection(set(geno.columns[1:])))

a_smp_i = np.array([smp2idx[k] for k,v in geno_at_marker.to_dict().items() if k in samples_shared and v == "B"])
b_smp_i = np.array([smp2idx[k] for k,v in geno_at_marker.to_dict().items() if k in samples_shared and v == "D"])

a_spectra_sum = np.sum(spectra[a_smp_i], axis=0)
b_spectra_sum = np.sum(spectra[b_smp_i], axis=0)

mutation_comparison(b_spectra_sum, a_spectra_sum, mut2idx)