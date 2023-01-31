from cyvcf2 import VCF 
import pandas as pd 
import numpy as np 

PROJDIR = "/scratch/ucgd/lustre-work/quinlan/u1006375/proj-mutator-mapping"

vcf = VCF(f"{PROJDIR}/data/vcf/AllMouse.vcf_90_recalibrated_snps_raw_indels_reheader_PopSorted.PASS.vcf.gz", gts012=True)
smp2idx = dict(zip(vcf.samples, range(len(vcf.samples))))
idx2smp = {v:k for k,v in smp2idx.items()}

chrom, start, end = "chr6", 113328509, 113328510
slop = 250_000
region = f"{chrom}:{start - slop}-{end + slop}"

geno_header = ["marker"]
geno_header.extend(vcf.samples)

marker_header = ["marker", "chromosome", "Mb"]

geno_outfh = open(f"{PROJDIR}/data/genotypes/wild.geno", "w")
print (','.join(geno_header), file=geno_outfh)

marker_outfh = open(f"{PROJDIR}/data/genotypes/wild.markers", "w")
print (','.join(marker_header), file=marker_outfh)

for v in vcf(region):
    marker = f"{v.CHROM}_{v.start}"
    mb = v.start / 1_000_000
    gq = v.gt_quals 
    ad, rd = v.gt_alt_depths, v.gt_ref_depths
    td = ad + rd 
    ab = ad / td
    gts = v.gt_types 

    good_idxs = np.where((td >= 10) & (gq >= 20) & (gts != 3))[0]

    if not 0.1 < (np.sum(gts[good_idxs]) / good_idxs.shape[0]) < 0.9: continue

    v_info = [marker]
    for idx, smp in idx2smp.items():
        geno = gts[idx]
        v_info.append(geno)

    m_info = list(map(str, [marker, v.CHROM, mb]))
    v_info = list(map(str, v_info))
    
    print (','.join(v_info), file=geno_outfh)
    print (','.join(m_info), file=marker_outfh)
