import pandas as pd 
from mutyper import Ancestor

def annotate_with_kmer(row, ancestor):
    chrom = row["chromosome"]
    ref, alt = row["ref"], row["alt"]
    pos = row["start"]
    kmer = ancestor.mutation_type(chrom, pos, ref, alt)
    kmer_form = None
    if None in kmer:
        kmer_form = "NNN>NNN"
    else:
        kmer_form = ">".join(kmer)
    return kmer_form

PROJDIR = "/scratch/ucgd/lustre-work/quinlan/u1006375/proj-mutator-mapping"

ref_fh = f"{PROJDIR}/data/ref/hg19.masked.fa"
mut_fh = f"{PROJDIR}/data/mutations/ceph/second_gen.dnms.txt"
smp_fh = f"{PROJDIR}/data/dbgap/phs001872.v1.pht009365.v1.p1.CEPH_Utah_Sample.MULTI.txt.gz"

mutations = pd.read_csv(mut_fh, sep="\t")
ancestor = Ancestor(ref_fh, sequence_always_upper=True, k=3)
sample_mapping = pd.read_csv(smp_fh, sep="\t", skiprows=10)

smpdict = dict(zip(sample_mapping["SUBJECT_ID"], sample_mapping["SAMPLE_ID"]))

mutations["chromosome"] = mutations["chrom"].apply(lambda c: f"chr{c}")
mutations = mutations[mutations["mut"] != "indel"]

mutations["kmer"] = mutations.apply(lambda row: annotate_with_kmer(row, ancestor), axis=1)
mutations["orig_sample_id"] = mutations["new_sample_id"].apply(lambda s: smpdict[s])

new_df = []
for i,row in mutations.iterrows():
    sample, haplotype = row["orig_sample_id"], row["phase"]
    if haplotype == "na": continue
    haplotype = "P" if haplotype == "paternal" else "M"
    true_sample = f"{sample}_{haplotype}"
    new_df.append({
        "sample": true_sample,
        "kmer": row["kmer"],
        "count": 1,
    })

new_df = pd.DataFrame(new_df)
print (new_df)
