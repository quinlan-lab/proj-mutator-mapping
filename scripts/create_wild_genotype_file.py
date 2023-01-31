from cyvcf2 import VCF 
import pandas as pd 
import numpy as np 

PROJDIR = "/scratch/ucgd/lustre-work/quinlan/u1006375/proj-mutator-mapping"

vcf = VCF(f"{PROJDIR}/data/vcf/wild.chr6.100Mb.120Mb.phased.vcf.gz", gts012=True)
smp2idx = dict(zip(vcf.samples, range(len(vcf.samples))))
idx2smp = {v:k for k,v in smp2idx.items()}

chrom, start, end = "chr6", 113328509, 113328510
slop = 250_000
region = f"{chrom}:{start - slop}-{end + slop}"

geno_header = ["marker"]
for s in vcf.samples:
    geno_header.append(f"{s}_A")
    geno_header.append(f"{s}_B")

marker_header = ["marker", "chromosome", "Mb"]

geno_outfh = open(f"{PROJDIR}/data/genotypes/wild.geno", "w")
print (','.join(geno_header), file=geno_outfh)

marker_outfh = open(f"{PROJDIR}/data/genotypes/wild.markers", "w")
print (','.join(marker_header), file=marker_outfh)

res = []

for v in vcf(region):
    marker = f"{v.CHROM}_{v.start}"
    mb = v.start / 1_000_000
    gts = np.array(v.genotypes)[:, :-1] 

    v_info = [marker]
    for idx, smp in idx2smp.items():
        hap_a_geno, hap_b_geno = gts[idx]
        for hg in (hap_a_geno, hap_b_geno):
            v_info.append("A" if hg == 0 else "B")
        if v.start == 113328509 and v.end == 113328510:
            res.append({
                'sample': smp,
                'genotype': sum(gts[idx]),
            })

    m_info = list(map(str, [marker, v.CHROM, mb]))
    v_info = list(map(str, v_info))
    
    print (','.join(v_info), file=geno_outfh)
    print (','.join(m_info), file=marker_outfh)

res_df = pd.DataFrame(res)
res_df['Subspecies'] = res_df['sample'].apply(lambda s: s.split('_')[0])
res_df_grouped = res_df.groupby('Subspecies').agg({'genotype': lambda g: sum(g) / (len(g) * 2)}).reset_index()
print (res_df_grouped)
