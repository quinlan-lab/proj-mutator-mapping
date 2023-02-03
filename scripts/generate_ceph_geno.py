# simple method to identify sites in the F1 generation that
# are present on either the paternal or maternal haplotype inherited
# by the F1 individual. these are sites that we can then use to
# test for potential mutator allele activity (or linkage to mutator
# alleles). the basic strategy is as follows:
# find alleles that are either HET or HOM_ALT in one parent and
# HOM_REF in the other, and inherited as HET by the F1. by definition
# (in the absence of genotype error) these variants are on the parental
# haplotype that donated them. we could also do something more nuanced
# by looking for sites that are HET in one parent and HET in the
# other, and identifying which parent donated the ALT allele using linkage,
# but the first, simpler strategy should produce enough sites that we can
# move forward.

from cyvcf2 import VCF
from peddy import Ped
import numpy as np
import pandas as pd
import tqdm
import argparse

PROJDIR = "/scratch/ucgd/lustre-work/quinlan/u1006375/proj-mutator-mapping"


def main(args):

    vcf_fh = f"/scratch/ucgd/lustre/UCGD_Datahub/Repository/AnalysisData/2018/A523/18-01-23_WashU-Yandell-CEPH_Sent/UCGD/GRCh37/VCF/Complete/18-01-23_WashU-Yandell-CEPH_Sent_Final_1529614168.vcf.gz"
    ped_fh = f"{PROJDIR}/data/ped/ceph.ped"

    # read in mapping from dbgap to true IDs
    smp_fh = f"{PROJDIR}/data/dbgap/phs001872.v1.pht009365.v1.p1.CEPH_Utah_Sample.MULTI.txt.gz"
    sample_mapping = pd.read_csv(smp_fh, sep="\t", skiprows=10)

    # create mapping from dbgap to true
    smpdict = dict(zip(sample_mapping["SUBJECT_ID"], sample_mapping["SAMPLE_ID"]))
    smpdict = {str(k):v for k,v in smpdict.items()}

    vcf = VCF(vcf_fh, gts012=True)
    ped = Ped(ped_fh)

    # define various groups of samples
    samples = [s for s in ped.samples()]
    kids = [k for k in samples if k.mom is not None and k.dad is not None]
    p0 = [k for k in samples if k.mom is None and k.dad is None]
    f1 = [
        k for k in kids
        if all([c is None for c in (
            k.mom.mom,
            k.mom.dad,
            k.dad.mom,
            k.dad.dad,
        )])
    ]
    f2 = [k for k in kids if k not in f1]

    p0_formatted, f1_formatted, f2_formatted = [], [], []
    for s in samples:
        if s.sample_id in smpdict:
            if s in p0:
                p0_formatted.append(str(smpdict[s.sample_id]))
            elif s in f1:
                f1_formatted.append(str(smpdict[s.sample_id]))
            elif s in f2:
                f2_formatted.append(str(smpdict[s.sample_id]))
            else: continue

    vcf.set_samples(f1_formatted + f2_formatted)

    chroms, chromlens = vcf.seqnames, vcf.seqlens
    chrom2len = dict(zip(chroms, chromlens))

    smp2idx = dict(zip(vcf.samples, range(len(vcf.samples))))

    p0_idxs = np.array([smp2idx[s] for s in p0_formatted if s in smp2idx])
    f1_idxs = np.array([smp2idx[s] for s in f1_formatted if s in smp2idx])
    f2_idxs = np.array([smp2idx[s] for s in f2_formatted if s in smp2idx])

    res = []

    window_size = 100_000

    geno_header = ["marker"]
    uniq_haps = []
    for s in f2:
        if s.sample_id not in smpdict: continue
        for h in ("P", "M"):
            uniq_haps.append(f"{smpdict[s.sample_id]}_{h}")
    geno_header.extend(uniq_haps)

    marker_header = ["marker", "chromosome", "Mb"]

    geno_outfh = open(f"{PROJDIR}/data/genotypes/per-chrom/{args.chrom}.ceph.geno", "w")
    marker_outfh = open(f"{PROJDIR}/data/genotypes/per-chrom/{args.chrom}.ceph.markers", "w")

    print (",".join(geno_header), file=geno_outfh)
    print (",".join(marker_header), file=marker_outfh)
    
    chrom_length = chrom2len[args.chrom]

    windows = np.arange(1, chrom_length, window_size)
    window_s, window_e = windows[:-1], windows[1:]

    for s, e in tqdm.tqdm(zip(window_s, window_e)):
        counted_in_window = 0
        for v in vcf(f"{args.chrom}:{s}-{e}"):

            if counted_in_window > 0: break

            if v.var_type != "snp": continue
            if v.call_rate != 1: continue
            if len(v.ALT) > 1: continue

            marker = f"{v.CHROM}_{v.POS}"

            marker_vals = [
                marker,
                v.CHROM,
                str(int(v.POS) / 1_000_000),
            ]
            print (",".join(marker_vals), file=marker_outfh)

            gts = v.gt_types
            if np.sum(gts[f1_idxs]) == 0: continue

            geno_vals = ["_".join([v.CHROM, str(v.POS)])]

            # by default, assume that both haplotypes have the REF
            # allele
            p_genotype, m_genotype = "A", "A"
            for s in f2:
                if s.sample_id not in smpdict:
                    continue
                true_s = str(smpdict[s.sample_id])
                
                kid_i = smp2idx[true_s]
                kid_gt = gts[kid_i]
                true_d, true_m = (
                    str(smpdict[s.dad.sample_id]),
                    str(smpdict[s.mom.sample_id]),
                )
                if true_d not in smp2idx or true_m not in smp2idx:
                    geno_vals.extend(["U", "U"])
                    continue
                dad_i, mom_i = smp2idx[true_d], smp2idx[true_m]
                dad_gt, mom_gt = gts[dad_i], gts[mom_i]
                if dad_gt > 0 and mom_gt > 0:
                    geno_vals.extend(["U", "U"])
                    continue

                # if the dad has an ALT at this site and the kid has
                # the ALT as well, then the paternal haplotype carries
                # the ALT and the maternal haplotype carries the REF
                if dad_gt > 0:# and kid_gt > 0:
                    p_genotype = "B"
                # if mom has the ALT and kid also has the ALT...
                elif mom_gt > 0:# and kid_gt > 0:
                    m_genotype = "B"
                # then, output this haplotype's genotypes at the site
                geno_vals.extend([p_genotype, m_genotype])
            if len(geno_vals) != len(geno_header): print ("ERROR")
            print (','.join(geno_vals), file=geno_outfh)
            counted_in_window += 1

if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--chrom")
    args = p.parse_args()
    main(args)