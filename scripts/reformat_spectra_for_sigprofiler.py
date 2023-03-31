import pandas as pd
import argparse


def revcomp(kmer: str) -> str:
    rc_dict = {"A":"T", "T":"A", "C":"G", "G":"C"}
    new_kmer = ""
    for k in list(kmer)[::-1]:
        new_kmer += rc_dict[k]
    return new_kmer

def find_central_mut(kmer: str) -> str:
    orig, new = kmer.split('>')
    fp, tp = orig[0], orig[2]
    central_orig = orig[1]
    central_new = new[1]
    return ">".join([central_orig, central_new])


def convert_kmer_to_sig_profiler_format(k):
    """
    convert the mutation types from in the BXD data
    to a format acceptable to SigProfilerExtractor. namely,
    convert A>T mutations to their reverse complements and
    format mutations as a A[B>D]C format, where A is the five-prime
    nucleotide, B>D is the mutation, and C is the three-prime nucleotide.
    """
    base_mut = find_central_mut(k)
    fprime, tprime = k[0], k[2]
    if base_mut[0] == "A":
        base_mut = "{}>{}".format(revcomp(base_mut.split('>')[0]),
                                  revcomp(base_mut.split('>')[1]))
        fprime, tprime = revcomp(k[0]), revcomp(k[2])
    new_mutation_type = "{}[{}]{}".format(fprime, base_mut, tprime)
    return new_mutation_type

def convert_kmer_to_sig_profiler_format_b(k):
    """
    convert the mutation types from in the BXD data
    to a format acceptable to SigProfilerExtractor. namely,
    convert A>T mutations to their reverse complements and
    format mutations as a A[B>D]C format, where A is the five-prime
    nucleotide, B>D is the mutation, and C is the three-prime nucleotide.
    """
    base_mut = find_central_mut(k)
    trinucleotide = k.split(">")[0]
    
    return base_mut, trinucleotide


def main(args):

    mutations = pd.read_csv(args.mutations)

    group_cols = [
        "sample",
        'kmer',
    ]

    # convert to wide-form dataframe
    mutations_tidy = mutations.groupby(group_cols).count().reset_index()
    mutations_tidy = mutations_tidy[group_cols + ["count"]]
    #mutations_tidy['Mutation type'] = mutations_tidy['kmer'].apply(lambda k: find_central_mut(k))
    #mutations_tidy['Trinucleotide'] = mutations_tidy['kmer'].apply(lambda k: k.split(">")[0])
    mutations_tidy['Mutation Types'] = mutations_tidy['kmer'].apply(lambda k: convert_kmer_to_sig_profiler_format(k))
    mutations_tidy.drop(columns=["kmer"], inplace=True)
    # # make a new column in the dataframe with the reformatted
    # # SigProfilerExtractor mutation type
    # mutations_tidy['Mutation Types'] = mutations_tidy['kmer'].apply(
    #     lambda k: convert_kmer_to_sig_profiler_format(k))

    groupcols = ["Mutation Types"]
    # convert to sigprofiler format
    mutations_sigpro = mutations_tidy.pivot(index=groupcols, columns="sample").reset_index()
    columns = groupcols.copy()
    columns.extend([c[1] for c in mutations_sigpro.columns if c[0] not in columns])
    mutations_sigpro.columns = columns 

    mutations_sigpro.fillna(value=0).to_csv(args.out, sep="\t", index=False)

if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument(
        "--mutations",
        required=True,
        help="""Path to CSV file containing per-mutation data.""",
    )
    p.add_argument(
        "--out",
        required=True,
    )
    args = p.parse_args()
    main(args)