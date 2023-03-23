import csv
import argparse
from collections import defaultdict, Counter
import gzip 
import csv 
from bx.intervals.intersection import IntervalTree
from pyfaidx import Fasta
import numpy as np
import pandas as pd
import tqdm

def make_interval_tree(
    path: str,
    datacol: bool = False,
    delim: str = '\t',
) -> defaultdict(IntervalTree):
    """
    generate an interval tree from a simple BED
    file with format chr:start:end
    """

    tree = defaultdict(IntervalTree)
    f = gzip.open(path, 'rt') if path.endswith('.gz') else open(path, 'r')
    fh = csv.reader(f, delimiter=delim)
    print ("Constructing exclude tree.")
    for line in fh:
        if datacol:
            tree[line[0]].insert(int(line[1]), int(line[2]), {'other': line[3]})
        else:
            tree[line[0]].insert(int(line[1]), int(line[2]))

    return tree

def main(args):

    res = []

    exclude = make_interval_tree(args.exclude)

    ref = Fasta(args.ref)

    chroms = list(map(str, range(1, 20)))
    chroms = ['chr' + c for c in chroms]

    nucs = ["A", "C"]
    nuc2idx = dict(zip(nucs, range(len(nucs))))

    # read in file with coverage information from mosdepth
    sample = args.coverage.split('/')[-1].split('.')[0]

    genome_wide_composition = np.zeros(2)

    with open(args.coverage, "r") as f:
        csvf = csv.reader(f, delimiter='\t')
        for l in tqdm.tqdm(csvf):
            # skip header line
            if l[0] == "#chrom": continue
            # make sure we only look at autosomes
            if l[0] not in chroms: continue
            # on each line, get the interval and the count
            # of bases in the interval that were covered by 
            # at least 10 reads
            chrom, start, end, _, _, _, tenx, _ = l

            callable_kmers = np.zeros(2)

            # figure out the nucleotide composition of the entire
            # interval and increment the global array
            sequence = ref[chrom][int(start):int(end)].seq
            nuc_counts = Counter(sequence.upper())
            for n, c in nuc_counts.items():
                if n == "N": continue 
                # reverse complement if necessary
                if n == "T": n = "A"
                if n == "G": n = "C"
                ni = nuc2idx[n]
                callable_kmers[ni] += c

            # figure out how much of this interval overlaps the masked regions
            exclude_overlaps = exclude[chrom].find(int(start), int(end))
            for e in exclude_overlaps:
                if e is None: continue
                sequence = ref[chrom][e.start:e.end].seq
                nuc_counts = Counter(sequence.upper())
                # subtract counts of nucleotides in excluded regions
                for n, c in nuc_counts.items():
                    if n == "N": continue 
                    # reverse complement if necessary
                    if n == "T": n = "A"
                    if n == "G": n = "C"
                    ni = nuc2idx[n]
                    callable_kmers[ni] -= c

            # we'll assume that if X% of bases were covered by at least 
            # 10 reads, then we can simply multiply the global kmer composition
            # counts by X%
            tenx_frac = int(tenx) / (int(end) - int(start))
            callable_kmers *= tenx_frac

            genome_wide_composition += callable_kmers
        
        for n, ni in nuc2idx.items():
            res.append({
                "sample": sample,
                "nucleotide": n,
                "count": int(genome_wide_composition[ni]),
            })

    res_df = pd.DataFrame(res)
    res_df.to_csv(args.out, index=False)

if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument(
        '--coverage',
        required=True,
        help="""path to sample threshold BED files""",
    )
    p.add_argument("--ref")
    p.add_argument(
        '--exclude',
        required=True,
        help="""path to file containing regions we want to mask""",
    )
    p.add_argument(
        '--out',
        required=True,
        help="""path to output file""",
    )
    args = p.parse_args()
    main(args)
