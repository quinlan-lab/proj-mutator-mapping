import pandas as pd
import numpy as np
import sys
import argparse
from skbio.stats.composition import clr

sys.path.insert(0, '/Users/tomsasani/quinlanlab/proj-mutator-mapping/')
from ihd.utils import compute_spectra

def main(args):

    # read in callable kmer content in each strain
    kmer_content = pd.read_csv(args.kmer_content)
    # read in sample genotypes at each marker
    genos = pd.read_csv(args.genos).set_index("marker")
    # define markers to use for the chr4 and chr6 mutator loci
    markers = np.array(["rs27509845", "rs46276051"])
    # map sample names to genotypes at each marker
    geno_at_marker = genos.loc[markers].to_dict()

    # read in annotated mutation data
    mutations = pd.read_csv(args.mutations)

    # compute mutation spectra in each sample using the specified k-mer context
    samples, mutation_types, spectra = compute_spectra(mutations, k=args.k)

    # map samples and mutation types to their indices in the mutation spectrum array
    smp2idx = dict(zip(samples, range(len(samples))))
    mut2idx = dict(zip(mutation_types, range(len(mutation_types))))
    samples_shared = list(set(samples).intersection(set(genos.columns)))

    # map samples to relevant metadata
    smp2generations = dict(zip(mutations['sample'], mutations['n_generations']))
    smp2epoch = dict(zip(mutations['sample'], mutations['true_epoch']))

    # map samples to their genotypes at both of the mutator loci markers
    smp2genotype = {}
    for s in samples_shared:
        sorted_markers = []
        for m in markers:
            sorted_markers.append(geno_at_marker[s][m])
        marker_genos = "-".join(sorted_markers)
        smp2genotype[s] = marker_genos

    # convert counts of each mutation type to fractions
    spectra_fracs = spectra / np.sum(spectra, axis=1)[:, np.newaxis]
    spectra_clr = clr(spectra_fracs)

    tidy_df = []
    # loop over strains
    for s in samples_shared:
        si = smp2idx[s]
        smp_geno = smp2genotype[s]
        smp_gens = smp2generations[s]
        if smp_gens < 0: continue
        # only consider samples with known homozygous genotypes at each site
        if "H" in smp_geno: continue
        # get counts of callable nucleotides in the strain
        callable_kmers = kmer_content[kmer_content["GeneNetwork name"] == s]

        for m, mi in mut2idx.items():
            # for each mutation type, get the "base" nucleotide
            base_nuc = m.split(">")[0]
            if base_nuc == "CpG": base_nuc = "C"
            # figure out the number of callable "base" nucleotides in that strain
            n_callable_bp = 1
            if args.k == 1:
                n_callable_bp = callable_kmers[callable_kmers["nucleotide"] == base_nuc]["count"].values[0]
            # convert counts of mutations to a rate
            rate = spectra[si, mi] / smp2generations[s] / n_callable_bp
            # append a bunch of info to the tidy dataframe
            vals = {
                'sample': s,
                'mut': m,
                'Mutation type': m,
                'Count': spectra[si, mi],
                'Total': np.sum(spectra, axis=1)[si],
                'Fraction': spectra_fracs[si, mi],
                'CLR_fraction': spectra_clr[si, mi],
                'Epoch': smp2epoch[s],
                "Generations": smp2generations[s],
                "callable_nucleotides": n_callable_bp,
                'ADJ_AGE': n_callable_bp * smp2generations[s],
                'Rate': rate,
                "is_ca": 1 if m == "C>A" else 0,
                "Haplotypes": smp_geno,
                'Haplotype_A': smp_geno.split('-')[0],
                'Haplotype_B': smp_geno.split('-')[1],
            }
            tidy_df.append(vals)
                

    tidy_df = pd.DataFrame(tidy_df)
    tidy_df.to_csv(args.out, index=False)

if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument(
        "--kmer_content",
        help="""counts of callable base pairs in each strain, stratified by nucleotide context""",
    )
    p.add_argument(
        "--genos",
        help="genotypes in each strain at each marker",
    )
    p.add_argument(
        "--mutations",
        help="""annotated and filtered mutation data in each strain""",
    )
    p.add_argument(
        "--out",
        help="""name of output file with tidy mutation spectra""",
    )
    p.add_argument(
        "-k",
        help="""kmer context in which to compute mutation spectra""",
        default=1,
        type=int,
    )
    args = p.parse_args()
    main(args)