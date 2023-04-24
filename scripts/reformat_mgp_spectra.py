import argparse
import pandas as pd
import numpy as np
from typing import Dict, List
from cyvcf2 import VCF

REV_COMP = {
    "A": "T",
    "T": "A",
    "C": "G",
    "G": "C",
}


def convert_to_mutation(row) -> str:
    """Convert 3-mer mutations (as reported in the Dumont 2019
    supplementary information) in each MGP strain to a 1-mer
    mutation type and collapse by strand complement such that all
    mutations have a "base" of either C or A.

    Args:
        row (_type_): Pandas row object.

    Returns:
        mutation (str): 1-mer mutation type string.
    """
    orig, new = row["ancestral"], row["derived"]
    fp, tp = row["5prime_flank"], row["3prime_flank"]
    mutation = None
    # if we need to reverse complement
    if orig in ("G", "T"):
        orig, new = REV_COMP[orig], REV_COMP[new]
        fp, tp = REV_COMP[fp], REV_COMP[tp]
    # treat CpG>TpG separately
    if orig == "C" and new == "T" and tp == "G":
        mutation = r"$\rightarrow$".join(["CpG", "TpG"])
    else:
        mutation = r"$\rightarrow$".join([orig, new])
    if mutation is None: print (row["FocalStrain"], row["chr"], row["pos"])
    return mutation


def get_genotype_at_site(vcf: VCF, pos: str) -> pd.DataFrame:
    """Get genotypes for every sample in the provided VCF and add
    to a pandas dataframe. 

    Args:
        vcf (VCF): Sanger Mouse Genomes Project VCF file
        pos (str): Position in GRCm38 coordinates

    Returns:
        pd.DataFrame: Dataframe with three columns (strain, pos, genotype)
    """
    genos = []
    for v in vcf(pos):
        # access ALT and REF depths from DP4 format field
        dp4 = v.format("DP4")
        ad, rd = np.sum(dp4[:, :2], axis=1), np.sum(dp4[:, 2:], axis=1)
        td = ad + rd
        gts = v.gt_types
        # loop over samples and collate genotypes
        for si, s in enumerate(vcf.samples):
            gt = -1
            # require total depth in sample to be >= 10
            if td[si] >= 10:
                gt = gts[si]
            genos.append({
                "strain": s,
                "genotype": "B" if gt == 0 else "D" if gt == 2 else "U"
            })

    # remove samples with unknown/low quality genotypes at the site and 
    # create a dataframe mapping samples to their genotypes at the locus
    geno_df = pd.DataFrame(genos).query('genotype != "U"')
    geno_df["pos"] = pos

    return geno_df


def combine_genotypes(genotypes_at_pos: List[pd.DataFrame], n_sites: int = 5) -> Dict[str, str]:
    """Combine a list of dataframes containing sample genotypes at `n_sites` total
    sites and filter to the samples that are homogenous for a single genotype at
    all `n_sites`. Then, create a dictionary mapping those samples to their
    homogenous genotypes.

    Args:
        genotypes_at_pos (List[pd.DataFrame]): List of dataframes containing mappings \
            of samples to genotypes at `n_sites` total sites.

        n_sites (int): Number of sites at which samples were genotyped.

    Returns:
        Dict[str, str]: Dictionary mapping samples to genotypes, only including those \
            samples that have the same genotype at every site.
    """
    
    geno_df = pd.concat(genotypes_at_pos)
    geno_df_agg = geno_df.groupby("strain").agg({"genotype": lambda g: homogenize_genotypes(g, n_sites=n_sites)}).reset_index()
    smp2geno = dict(zip(geno_df_agg["strain"], geno_df_agg["genotype"]))
    return smp2geno


def homogenize_genotypes(genotypes: List[str], n_sites: int = 5) -> str:
    """Combine sample genotypes across `n_sites` such that we only
    consider samples to be D- or B-like if they possess D or B alleles
    at every one of the `n_sites`. Otherwise, they're considered "heterozygous."

    Args:
        genotypes (List[str]): List of genotypes at `n_sites` sites.

        n_sites (int, optional): Number of sites at which we expect samples to have genotypes. Defaults to 5.

    Returns:
        str: Consensus genotype.
    """
    
    if len(genotypes) < n_sites: return "U"
    if len(set(genotypes)) > 1: return "H"
    if all([g == "D" for g in genotypes]): return "D"
    elif all([g == "B" for g in genotypes]): return "B"
    else: return "U"


def main(args):

    # read in MGP VCF
    vcf = VCF(args.vcf, gts012=True)

    # figure out which MGP strains possess D- or B-like alleles
    # at the markers on chromosome 4 and chromosome 6
    chr4_marker_pos = "4:116750874-116750875"
    chr6_marker_pos = "6:114052413-114052414"

    ### actual missense positions
    chr4_marker_pos = [
        "4:116814337-116814338",
        "4:116814393-116814394",
        "4:116815657-116815658",
        "4:116817415-116817416",
        "4:116817418-116817419",
    ]
    chr6_marker_pos = [
        "6:113328509-113328510",
    ]

    # since there are 5 possible missense mutations at the chromosome 4
    # mutator locus, we'll first identify the MGP strains that are D-like and
    # B-like at all of the sites.
    chr4_genos_df = []
    for chr4_pos in chr4_marker_pos:
        chr4_genos_df.append(get_genotype_at_site(vcf, chr4_pos))
    chr4_genos = combine_genotypes(chr4_genos_df, n_sites = len(chr4_marker_pos))

    chr6_genos_df = []
    for chr6_pos in chr6_marker_pos:
        chr6_genos_df.append(get_genotype_at_site(vcf, chr6_pos))
    chr6_genos = combine_genotypes(chr6_genos_df, n_sites = len(chr6_marker_pos))

    # read in numbers of callable base pairs in each strain
    dumont_callable = pd.read_excel(
        args.dumont_xls,
        sheet_name="TableS3",
        header=2,
    )

    cols2use = ["Strain"]
    cols2use.extend(["nC", "nT", "nA", "nG"])

    # reformat column names
    dumont_callable = dumont_callable[cols2use]
    dumont_callable['strain'] = dumont_callable['Strain'].apply(
        lambda x: x.replace('/', '_'))

    # read in per-sample mutation data in MGP strains
    dumont_long = pd.read_excel(args.dumont_xls, sheet_name="TableS4", header=2)
    # add formatted mutation column to each row
    dumont_long["Mutation type"] = dumont_long.apply(lambda row: convert_to_mutation(row), axis=1)
    # figure out total counts of each mutation type in each MGP strain
    dumont_long_grouped = dumont_long.groupby([
        "FocalStrain",
        "Mutation type",
    ]).size().reset_index().rename(columns={
        "FocalStrain": "strain",
        0: "Count",
    })

    # merge per-sample mutation counts with callable base pairs
    dumont_merged = dumont_long_grouped.merge(dumont_callable, on="strain")
    # annotate merged dataframe with information about haplotypes at each marker
    dumont_merged["Haplotypes"] = dumont_merged["strain"].apply(lambda s: "-".join([chr4_genos[s], chr6_genos[s]]))
    # figure out total numbers of callable cytosines and adenines in each strain
    dumont_merged["CALLABLE_C"] = dumont_merged["nC"] + dumont_merged["nG"]
    dumont_merged["CALLABLE_A"] = dumont_merged["nA"] + dumont_merged["nT"]
    # and drop superfluous columns
    dumont_merged.drop(columns=["nC", "nT", "nA", "nG"], inplace=True)

    # convert mutation counts into rates, taking into account the nucleotide context
    # in which each 1-mer occurred
    dumont_merged["Rate"] = dumont_merged.apply(
        lambda row: row["Count"] / row["CALLABLE_C"] if row["Mutation type"].
        startswith("C") else row["Count"] / row["CALLABLE_A"],
        axis=1,
    )
    # sum the context-adjusted rates in each strain
    total_rates = dumont_merged.groupby(["strain"]).agg({"Rate": sum}).reset_index().rename(columns={"Rate": "Total rate"})

    # create a tidy dataframe with 1-mer and total rates in each strain
    dumont_tidy = dumont_merged.merge(total_rates, on="strain")
    # convert rates to fractions as well
    dumont_tidy["Fraction"] = dumont_tidy["Rate"] / dumont_tidy["Total rate"]
    # filter the tidy dataframe to only include strains with homogenously D or B genotypes
    # at both of the mutator loci
    dumont_tidy = dumont_tidy[dumont_tidy["Haplotypes"].isin(["B-B", "B-D", "D-B", "D-D"])]
    dumont_tidy["is_ca"] = dumont_tidy["Mutation type"].apply(lambda m: 1 if m == r"C$\rightarrow$A" else 0)
    dumont_tidy["Haplotype_A"] = dumont_tidy["Haplotypes"].apply(lambda h: h.split('-')[0])
    dumont_tidy["Haplotype_B"] = dumont_tidy["Haplotypes"].apply(lambda h: h.split('-')[1])

    dumont_tidy.to_csv(args.out, index=False)


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument(
        "--dumont_xls",
        help="""strain-private substitution data from Dumont et al. (2019) \
                            as it appears in Table 3 of the Supplementary Data""",
    )
    p.add_argument(
        "--vcf",
        help="VCF file containing relevant mutation data in MGP strains",
    )
    p.add_argument(
        "--out",
        help="""name of output file with MGP spectra""",
    )
    args = p.parse_args()
    main(args)