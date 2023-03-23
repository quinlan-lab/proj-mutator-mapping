import pandas as pd 

sv_fh = "ftp://ftp.ebi.ac.uk/pub/databases/mousegenomes//REL-1606-SV/REL-1606-SV/mgpv5.SV_deletions.bed.gz"
sv_df = pd.read_csv(sv_fh, sep="\t", dtype={"#chr": str, "start": int, "end": int})


sv_df = sv_df[(sv_df["#chr"] == "6") & ( (sv_df["start"] >= 109_000_000) & (sv_df["end"] <= 119_000_000) )]

columns = ["chrom", "start", "end", "C57BL_6NJ", "DBA_2J"]

unk = ".;.;.;."


sv_df = sv_df[(sv_df["C57BL_6NJ"] != unk) | (sv_df["DBA_2J"] != unk)]
sv_df["chrom"] = sv_df["#chr"].apply(lambda c: f"chr{c}")

sv_df[columns].to_csv("dels.tsv", sep="\t", index=False, header=False)

refseq_fh = "/Users/tomsasani/Downloads/mm10.refseq.tsv"
refseq_df = pd.read_csv(refseq_fh, sep="\t")
refseq_df = refseq_df[(refseq_df["chrom"] == "chr6") & (refseq_df["txStart"] >= 109_000_000) & (refseq_df["txEnd"] <= 119_000_000)]

new_refseq_df = []
for i,row in refseq_df.iterrows():
    exon_starts, exon_ends = row["exonStarts"].split(","), row["exonEnds"].split(",")
    for i, (s, e) in enumerate(zip(exon_starts, exon_ends)):
        new_refseq_df.append({
            "chrom": row["chrom"],
            "start": s,
            "end": e,
            "exon_num": i + 1,
            "gene": row["name2"],
            "transcript_name": row["name"],

        })

new_refseq_df = pd.DataFrame(new_refseq_df)
new_refseq_df.to_csv("new_refseq.tsv", sep="\t", index=False, header=False)
