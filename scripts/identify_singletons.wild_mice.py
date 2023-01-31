import numpy as np
import argparse
from collections import defaultdict
from cyvcf2 import VCF
from mutyper.ancestor import Ancestor
import doctest
import pandas as pd
from bx.intervals.intersection import Interval, IntervalTree
import gzip 
import csv 
from typing import Tuple

def get_good_idxs(
    gts: np.ndarray,
    gq: np.ndarray,
    td: np.ndarray,
    min_dp: int = 10,
    min_gq: int = 20,
) -> np.ndarray:
    """
    get a list of indices corresponding to samples
    that meet a set of reasonable filters.

    gts: 	array of sample genotypes
    gq: 	array of sample genotype qualities 
    td: 	array of sample depths

    """

    UNK, HOM_REF, HET, HOM_ALT = range(-1, 3)

    # get indices of non-unk genotypes
    called_gts = np.where(gts != UNK)[0]
    # and indices of genotypes that meet GQ and DP filters
    good_gq = np.where(gq >= min_gq)
    good_depth = np.where(td >= min_dp)

    # intersect all of the above indices together to get
    # the output indices that meet all filters
    good_sites = np.intersect1d(good_gq, good_depth)
    good_sites_not_unk = np.intersect1d(good_sites, called_gts)

    return good_sites_not_unk

def normalize_var(ref: str, alt: str) -> Tuple[str, str]:
    """
	method for normalizing a variant by hand
	"""

    ref_a = np.array(list(ref))
    alt_a = np.array(list(alt))

    diff_nucs = ref_a != alt_a
    diff_idxs = np.where(diff_nucs)[0]

    if diff_idxs.shape[0] == 0: return ('N', 'N')
    elif diff_idxs.shape[0] > 1: return ('N', 'N')
    else:
        return (ref_a[diff_idxs[0]], alt_a[diff_idxs[0]])


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


def main():

    # -----------
    # read in VCF
    # -----------
    PROJDIR = "/scratch/ucgd/lustre-work/quinlan/u1006375/proj-mutator-mapping"
    vcf = VCF(f"{PROJDIR}/data/vcf/AllMouse.vcf_90_recalibrated_snps_raw_indels_reheader_PopSorted.PASS.vcf.gz", gts012=True)

    # ------------
    # define some global variables that will come in handy
    # ------------
    UNK, HOM_REF, HET, HOM_ALT = range(-1, 3)

    ref_fh = "/scratch/ucgd/lustre-work/quinlan/u1006375/proj-mouse-spectra/data/ref/mm10.fa"
    ancestor = Ancestor(ref_fh, k=3, sequence_always_upper=True)

    # map sample indices to sample IDs and vice versa
    smp2idx = dict(zip(vcf.samples, range(len(vcf.samples))))
    idx2smp = {v:k for k,v in smp2idx.items()}

    # map "ancestries" (i.e., mouse subspecies) to sample indexes
    anc2idx = defaultdict(list)
    anc = [s.split('_')[0] for s in vcf.samples]
    smp2anc = dict(zip(vcf.samples, anc))
    for smp in smp2idx:
        anc = smp2anc[smp]
        idx = smp2idx[smp]
        anc2idx[anc].append(idx)

    # ------------
    # iterate over the VCF
    # ------------
    chrom, start, end = "chr6", 113328509, 113328510
    slop = 250_000
    region = f"{chrom}:{start - slop}-{end + slop}"
    vcf_h = vcf(region)

    # generate an exclude file if it's provided
    exclude_fh = f"/scratch/ucgd/lustre-work/quinlan/u1006375/proj-mouse-spectra/data/mm10.seg_dups.simple_repeats.merged.bed.gz"
    exclude = make_interval_tree(exclude_fh)

    res = []

    for i,v in enumerate(vcf_h):
        # check if the current variant is in the exclude file
        if len(exclude[v.CHROM].find(v.start, v.end)) > 0:
            continue
        if v.var_type != "snp": continue
        if '*' in v.ALT: continue
        if len(v.REF) > 1: continue
        if len(v.ALT) > 1: continue
        # require all samples to have a callable genotype
        if v.call_rate < 1: continue
    

        for alt_idx, alt in enumerate(v.ALT):
            if len(alt) > 1: continue

            ref_allele, alt_allele = v.REF, alt

            # don't consider indels
            if len(v.REF) != len(alt): continue

            # if the REF allele is longer than one nucleotide, but
            # the REF and ALT alleles are the same length, we almost
            # certainly need to normalize to become an SNV
            if len(v.REF) > 1 and (len(v.REF) == len(alt)):
                ref_allele, alt_allele = normalize_var(v.REF, alt)
                if 'N' in (ref_allele, alt_allele): continue

            # define the alternate genotype at this site.
            # if this is a multiallelic site, the genotype (integer)
            # corresponding to the first allele is 1, the genotype
            # corresponding to the second allele is 2, and so on.
            alt_gt = alt_idx + 1

            # access sample genotypes, excluding boolean flag indicating
            # whether sample is phased
            gts = np.array(v.genotypes)[:, :-1]
            if np.sum(gts) > 1: continue

            # the index of the reference allele is always 0
            ref_idx = 0

            # access ref and alt depths
            rd = v.format('AD')[:, ref_idx]
            ad = v.format('AD')[:, alt_idx + 1]
            rd[rd < 0] = 0
            ad[ad < 0] = 0

            td = rd + ad
            ab = ad / td
            gq = v.gt_quals

            # get a list of sample indices that meet quality thresholds
            # e.g., non-UNK, DP>=10, GQ>=20
            for hap_i in [0, 1]:
                good_idxs = get_good_idxs(
                    gts[:, hap_i],
                    gq,
                    td,
                    min_dp=10,
                    min_gq=20,
                )

                if good_idxs.shape[0] == 0: continue

                # get the mutation context for this sites using mutyper if it's a SNP
                # otherwise, we'll just report the kmer is an "indel"
                kmer = ancestor.mutation_type(v.CHROM, v.start, ref_allele,
                                            alt_allele)
                if None in kmer: kmer = "NNN>NNN"
                else: kmer = '>'.join(list(kmer))
                if kmer == "NNN>NNN": continue

                smps_with_mut = np.where(gts[:, hap_i] == 1)
                good_smps_with_mut = np.intersect1d(smps_with_mut, good_idxs)

                for idx in good_smps_with_mut:
                    sample = idx2smp[idx]
                    haplotype = "A" if hap_i == 0 else "B"
                    d = {
                        'Strain': f"{sample}_{haplotype}",
                        'kmer': kmer,
                        'chrom': v.CHROM,
                        'start': v.start,
                        'end': v.end,
                    }
                    res.append(d)

    df = pd.DataFrame(res)
    df.to_csv(f"{PROJDIR}/data/singletons/wild/annotated_filtered_singletons.csv", index=False)

if __name__ == "__main__":
    main()
