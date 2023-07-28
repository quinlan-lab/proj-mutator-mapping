import pandas as pd
from cyvcf2 import VCF
from bx.intervals.intersection import Interval, IntervalTree
import argparse

def main(args):

    refseq = pd.read_csv(args.refseq, sep="\t")

    start = 111_270_000
    buff = 10_000_000

    chrom, start, end = "chr6", start - buff, start + buff

    refseq = refseq[(refseq["chrom"] == chrom) & (refseq["txStart"] > start)
                    & (refseq["txStart"] < end)]

    vcf = VCF(args.vcf)

    exons, full_tx = IntervalTree(), IntervalTree()

    for i, row in refseq.iterrows():
        # limit to curated refseq 
        #if not row["name"].startswith("NM"): continue
        exon_starts, exon_ends = row["exonStarts"].split(','), row["exonEnds"].split(',')
        tx_start, tx_end = int(row["txStart"]), int(row["txEnd"])
        tx_start -= 10_000
        tx_end += 10_000
        gene_name = row["name2"]
        for exon_start, exon_end in zip(exon_starts, exon_ends):
            if exon_start == "" and exon_end == "": continue
            exons.insert(int(exon_start), int(exon_end), {"gene": gene_name, "exon_start": exon_start, "exon_end": exon_end})
        full_tx.insert(tx_start, tx_end, {"gene": gene_name, "tx_start": tx_start, "tx_end": tx_end})

    smp2idx = dict(zip(vcf.samples, range(len(vcf.samples))))
    DBA, C57 = smp2idx["DBA"], smp2idx["6NJ"]

    out_df = []

    for v in vcf(f"{chrom}:{start}-{end}"):
        gts = v.gt_types
        if gts[DBA] + gts[C57] == 0: continue
        v_start, v_end = v.start, v.end
        if v.INFO.get("SVTYPE") == "DEL":
            v_end -= v.INFO.get("SVLEN")
        exon_overlaps = exons.find(v_start, v_end)
        tx_overlaps = full_tx.find(v_start, v_end)
        overlapping_genes = list(set([o["gene"] for o in tx_overlaps]))
        
        if len(tx_overlaps) == 0 and len(exon_overlaps) == 0: continue
        is_in_exon = False
        if len(exon_overlaps) > 0: is_in_exon = True
        vals = {
            "start": v_start,
            "end": v_end,
            "sv_type": v.INFO.get("SVTYPE"),
            "sv_len": v.INFO.get("SVLEN"),
            "in_exon": int(is_in_exon),
            "genes": "-".join(overlapping_genes),
            "c57_gt": gts[C57],
            "dba_gt": gts[DBA],
        }
        out_df.append(vals)

    out_df = pd.DataFrame(out_df)
    out_df.to_csv(args.out, index=False)

if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument(
        "--refseq",
        help=
        """Path to TSV file containing RefSeq transcript information, downloaded from the UCSC table browser.""",
    )
    p.add_argument(
        "--vcf",
        help="""Path to SV VCF file from Ferraj et al. (2023)""",
    )
    p.add_argument(
        "--out",
        help="""Path to output CSV file.""",
    )
    args = p.parse_args()
    main(args)
