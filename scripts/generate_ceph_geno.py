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

PROJDIR = "/scratch/ucgd/lustre-work/quinlan/u1006375/proj-mutator-mapping"

vcf_fh = f"{PROJDIR}/data/vcf/ceph.vcf.gz"
ped_fh = f"{PROJDIR}/data/ped/ceph.ped"

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

samples2use = p0 + f1
vcf.set_samples(samples2use)

smp2idx = dict(zip(vcf.samples, range(len(vcf.samples))))

p0_idxs = np.array([smp2idx[s.sample_id] for s in p0])
f1_idxs = np.array([smp2idx[s.sample_id] for s in f1])

res = []

for v in vcf("chr19"):
    if v.var_type != "snp": continue
    if v.call_rate != 1: continue
    if len(v.ALT) > 1: continue

    marker = f"{v.CHROM}_{v.start}"

    gts = v.gt_types
    if np.sum(gts[p0_idxs]) == 0: continue
    if np.sum(gts[f1_idxs]) == 0: continue


    for s in f1:
        kid_i = smp2idx[s.sample_id]
        kid_gt = gts[kid_i]
        dad_i, mom_i = smp2idx[s.dad.sample_id], smp2idx[s.mom.sample_id]
        dad_gt, mom_gt = gts[dad_i], gts[mom_i]
        if dad_gt > 0 and mom_gt > 0: continue
        # by default, assume that both haplotypes have the REF
        # allele
        p_genotype, m_genotype = "A", "A"
        # if the dad has an ALT at this site and the kid has
        # the ALT as well, then the paternal haplotype carries
        # the ALT and the maternal haplotype carries the REF
        if dad_gt > 0 and kid_gt > 0:
            p_genotype = "B"
        # if mom has the ALT and kid also has the ALT...
        elif mom_gt > 0 and kid_gt > 0:
            m_genotype = "B"
        # then, output this haplotype's genotypes at the site
        for hap, gen in zip(
            (f"{s.sample_id}_P", f"{s.sample_id}_M"),
            (p_genotype, m_genotype),
        ):
            res.append({
                "chromosome": v.CHROM,
                "position": v.POS,
                "haplotype": hap,
                "genotype": gen,
            })
