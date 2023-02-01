from cyvcf2 import VCF 
import numpy as np
from bx.intervals.intersection import IntervalTree
from collections import defaultdict 
import gzip
import csv 
import pandas as pd


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
        chrom, start, end, genes = line
        uniq_genes = list(set(genes.split(',')))
        tree[chrom].insert(int(start), int(end), {'gene': ','.join(uniq_genes)})
        
    return tree

PROJDIR = "/scratch/ucgd/lustre-work/quinlan/u1006375/proj-mutator-mapping"

vcf = VCF(f"{PROJDIR}/data/vcf/combined/merged.regenotyped.filtered.vcf.gz", gts012=True)
idx2smp = dict(zip(range(len(vcf.samples)), vcf.samples))

gene_tree = make_interval_tree(f"{PROJDIR}/data/GRCm39.dna_repair.sorted.merged.bed")

DBA = "sample_4512-JFI-0334_DBA_2J_three_lanes_phased_possorted_bam"
C57 = "sample_4512-JFI-0333_C57BL_6J_two_lanes_phased_possorted_bam"

res = []

for v in vcf:
    if v.FILTER not in (None, "PASS"): continue 
    # only consider DELs
    #if v.INFO.get("SVTYPE") != "DEL": continue 
    if v.INFO.get('SVLEN') is not None:
        if not (250 <= v.INFO.get("SVLEN") <= 10_000_000): continue 
    if v.INFO.get("SU") < 5: continue 

    supporting_ev = v.format("SU")
    if np.sum(supporting_ev > 0) > 1: continue
    # require at least one of paired-read or split-read evidence
    #if v.INFO.get("PE") == 0 and v.INFO.get("SR") == 0: continue

    start, end = v.start, v.INFO.get("END")

    # check overlap against DNA repair genes 
    gene_overlap = gene_tree[v.CHROM].find(start, end)
    if len(gene_overlap) == 0: continue 

    # get sample with singleton 
    singleton_idxs = np.where(v.format("SU") > 0)[0]
    if singleton_idxs.shape[0] == 0: continue 

    for idx in singleton_idxs:
        singleton_smp = idx2smp[idx].split('.')[0]
        for d in gene_overlap:
            res.append({
                'chrom': v.CHROM,
                'start': start,
                'end': end,
                'svtype': v.INFO.get('SVTYPE'),
                'sample': singleton_smp,
                'pe': v.format("PE")[idx][0],
                'sr': v.format("SR")[idx][0],
                'su': v.format("SU")[idx][0],
                'gene': d['gene'],
            })

res_df = pd.DataFrame(res)
res_df.to_csv(f"{PROJDIR}/hq_singleton_svs.csv", index=False)


    