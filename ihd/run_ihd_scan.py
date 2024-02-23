import pandas as pd
import numpy as np
import argparse
import json
import sys
from utils import (
    compute_spectra,
    perform_permutation_test,
    perform_ihd_scan,
    compute_manual_chisquare,
    compute_manual_cosine_distance,
    get_covariate_matrix,
    calculate_covariate_by_marker,
)
from schema import IHDResultSchema, MutationSchema
import numba


def filter_mutation_data(
    mutations: pd.DataFrame,
    geno: pd.DataFrame,
) -> pd.DataFrame:
    # get unique samples in mutation dataframe
    samples = mutations["sample"].unique()
    # get the overlap between those and the sample names in the genotype data
    samples_overlap = list(set(samples).intersection(set(geno.columns)))

    if len(samples_overlap) == 0:
        print(
            """Sorry, no samples in common between mutation data
        and genotype matrix. Please ensure sample names are identical."""
        )
        sys.exit()

    # then subset the genotype and mutation data to include only those samples
    cols2use = ["marker"]
    cols2use.extend(samples_overlap)
    geno = geno[cols2use]

    mutations_filtered = mutations[mutations["sample"].isin(samples_overlap)]
    return mutations_filtered


def main(args):
    # read in JSON file with file paths
    config_dict = None
    with open(args.config, "rb") as config:
        config_dict = json.load(config)

    # read in genotype info
    geno = pd.read_csv(config_dict["geno"])
    markers = pd.read_csv(config_dict["markers"])

    markers2use = markers[markers["chromosome"] != "X"]["marker"].unique()
    geno = geno[geno["marker"].isin(markers2use)]

    # read in singleton data and validate with pandera
    mutations = pd.read_csv(args.mutations, dtype={"sample": str})
    MutationSchema.validate(mutations)

    mutations_filtered = filter_mutation_data(mutations, geno)

    # get a list of samples and their corresponding mutation spectra
    samples, mutation_types, spectra = compute_spectra(
        mutations_filtered,
        k=args.k,
        cpg=True,
    )
    print(
        f"""Using {len(samples)} samples
          and {int(np.sum(spectra))} total mutations."""
    )

    strata = np.ones(len(samples))
    if args.stratify_column is not None:
        sample2strata = dict(
            zip(
                mutations_filtered["sample"],
                mutations_filtered[args.stratify_column],
            )
        )
        strata = np.array([sample2strata[s] for s in samples])

    callable_kmer_arr = None
    if args.callable_kmers and args.k == 1:
        callable_kmer_arr = np.zeros(
            (len(samples), len(mutation_types)), dtype=np.int64
        )
        callable_kmers = pd.read_csv(args.callable_kmers)
        # NOTE: need to check schema of df
        for si, s in enumerate(samples):
            for mi, m in enumerate(mutation_types):
                base_nuc = m.split(">")[0] if m != "CpG>TpG" else "C"
                callable_k = callable_kmers[
                    (callable_kmers["sample"] == s)
                    & (callable_kmers["nucleotide"] == base_nuc)
                ]
                callable_kmer_arr[si, mi] = callable_k["count"].values[0]

    # convert string genotypes to integers based on config definition
    replace_dict = config_dict["genotypes"]
    geno_asint = geno.replace(replace_dict).replace({1: np.nan})

    if args.adj_region is not None:
        chrom = args.adj_region.split(":")[0]
        start, end = list(map(float, args.adj_region.split(":")[1].split("-")))
        # find markers within this region
        markers_to_filter = markers[
            (markers["chromosome"] == chrom)
            & (markers["Mb"] >= start)
            & (markers["Mb"] <= end)
        ]["marker"].unique()
        marker_idxs = geno_asint[
            geno_asint["marker"].isin(markers_to_filter)
        ].index.values
        geno_asint = geno_asint.iloc[geno_asint.index.difference(marker_idxs)]

    # convert genotype values to a matrix
    geno_asint_filtered_matrix = geno_asint[samples].values
    # get an array of marker names at the filtered genotyped loci
    markers_filtered = geno_asint["marker"].values

    # compute similarity between allele frequencies at each marker
    # genotype_similarity = compute_genotype_similarity(geno_asint_filtered_matrix)
    genotype_similarity = np.ones(geno_asint_filtered_matrix.shape[0])
    distance_method = compute_manual_cosine_distance
    if args.distance_method == "chisquare":
        distance_method = compute_manual_chisquare

    covariate_cols = []
    covariate_matrix = get_covariate_matrix(
        mutations_filtered,
        samples,
        covariate_cols=covariate_cols,
    )

    covariate_ratios = calculate_covariate_by_marker(
        covariate_matrix,
        geno_asint_filtered_matrix,
    )

    # compute the maximum cosine distance between groups of
    # haplotypes at each site in the genotype matrix
    focal_dists = perform_ihd_scan(
        spectra,
        geno_asint_filtered_matrix,
        genotype_similarity,
        covariate_ratios,
        distance_method=distance_method,
        adjust_statistics=False,
    )

    res_df = pd.DataFrame(
        {
            "marker": markers_filtered,
            "Distance": focal_dists,
            "k": args.k,
        }
    )
    IHDResultSchema.validate(res_df)

    numba.set_num_threads(args.threads)

    # then do permutations
    null_distances = perform_permutation_test(
        spectra,
        geno_asint_filtered_matrix,
        genotype_similarity,
        covariate_ratios,
        strata,
        distance_method=distance_method,
        n_permutations=args.permutations,
        progress=args.progress,
        adjust_statistics=False,
    )

    # compute the Nth percentile of the maximum distance
    # distribution to figure out the distance thresholds
    for pctile in (20, 5, 1):
        score_pctile = np.percentile(null_distances, 100 - pctile, axis=0)
        res_df[f"{100 - pctile}th_percentile"] = score_pctile

    # for each chromosome, compute the specified confidence interval around
    # the peak observed distance
    geno_asint_filtered_merged = geno_asint.merge(markers, on="marker")

    # combined_conf_int_df = []

    # conf_int_chroms = ["4", "6"]
    # for chrom, chrom_df in geno_asint_filtered_merged.groupby("chromosome"):
    #     if chrom not in conf_int_chroms:
    #         continue

    #     chrom_genotype_matrix = chrom_df[samples].values
    #     # compute confidence intervals on the chromosome
    #     conf_int_lo, conf_int_hi, peak_markers = calculate_confint(
    #         spectra,
    #         chrom_genotype_matrix,
    #         covariate_matrix,
    #         distance_method=distance_method,
    #         adjust_statistics=True,
    #         conf_int=90.,
    #         n_permutations=1_000,
    #     )

    #     conf_int_df = chrom_df.iloc[np.array([conf_int_lo, conf_int_hi])]
    #     combined_conf_int_df.append(conf_int_df[["chromosome", "marker", "Mb"]])

    # combined_conf_int_df = pd.concat(combined_conf_int_df)

    # combined_conf_int_df.to_csv(f"{args.out}.ci.csv", index=False)

    res_df.to_csv(args.out, index=False)


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument(
        "--mutations",
        type=str,
        help="Path to mutation data in CSV format.",
    )
    p.add_argument(
        "--config",
        type=str,
        help="Path to config file in JSON format.",
    )
    p.add_argument(
        "--out",
        help="Path in which to store the results of the IHD scan.",
    )
    p.add_argument(
        "-k",
        type=int,
        default=1,
        help="k-mer context used to classify mutations. Default is 1.",
    )
    p.add_argument(
        "-permutations",
        type=int,
        default=1_000,
        help="Number of permutations to perform when calculating significance thresholds. Default is 1,000.",
    )
    p.add_argument(
        "-distance_method",
        default="cosine",
        type=str,
        help="""Method to use for calculating distance between aggregate spectra. Options are 'cosine' and 'chisquare', default is 'chisquare'.""",
    )
    p.add_argument(
        "-threads",
        default=1,
        type=int,
        help="""Number of threads to use during permutation testing step. Default is 1.""",
    )
    p.add_argument(
        "-progress",
        action="store_true",
        help="""Whether to output the progress of the permutation testing step (i.e., the number of completed permutations).""",
    )
    p.add_argument(
        "-callable_kmers",
        default=None,
        type=str,
        help="""Path to CSV file containing numbers of callable base pairs in each sample, stratified by nucleotide.""",
    )
    p.add_argument(
        "-stratify_column",
        default=None,
        type=str,
        help="""If specified, use this column to perform a stratified permutation test by only permuting BXDs within groups defined by the column to account for population structure.""",
    )
    p.add_argument(
        "-adj_region",
        default=None,
        type=str,
        help="""If specified, a chromosomal region (chr:start-end) that we should adjust for in our AMSD scans. Start and end coordinates should be specified in Mb. Default is None""",
    )
    args = p.parse_args()

    main(args)
