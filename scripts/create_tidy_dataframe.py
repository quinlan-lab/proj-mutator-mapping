import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import argparse
import json
import sys
from utils import compute_spectra

from schema import IHDResultSchema, MutationSchema

PROJDIR = "/Users/tomsasani/quinlanlab/proj-mutator-mapping"

def main():

    # read in singleton data and validate with pandera
    mutations = pd.read_csv(f"{PROJDIR}/data/mutations/bxd/annotated_filtered_singletons.csv")
    MutationSchema.validate(mutations)

    samples = mutations['sample'].unique()

    # for the null permutation test, shuffle the rows of the spectra
    # dataframe every time. otherwise keep it the same.
    samples, mutation_types, spectra = compute_spectra(mutations, k=1)

    spectra_sums = spectra

    spectra_df = pd.DataFrame(spectra, columns=mutation_types)
    spectra_df["sample"] = samples
    print (spectra_df)
    

if __name__ == "__main__":
    main()
