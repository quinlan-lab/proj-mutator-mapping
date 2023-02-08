import pandas as pd
import argparse

def main(args):
    res = []
    for fh in args.csvs:
        df = pd.read_csv(fh)
        res.append(df)
    res_df = pd.concat(res)

    metadata = pd.read_excel(args.metadata).dropna()

    res_df_merged = res_df.merge(
        metadata,
        left_on="sample",
        right_on="bam_name",
    )

    res_df_merged.to_csv(args.out, index=False)

if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument(
        "--csvs",
        nargs="+",
        help=
        "path to the CSVs containing callable kmer information \
            that you want to combine"                                     ,
    )
    p.add_argument(
        "--metadata",
        help="path to Excel file containing BXD metadata",
    )
    p.add_argument("--out", help="name of output file")
    args = p.parse_args()
    main(args)
