import pandas as pd 
import argparse

def main(args):
    combined = []
    for fh in args.genotypes:
        df = pd.read_csv(fh, skiprows=3)
        combined.append(df)

    combined = pd.concat(combined)

    combined.to_csv(args.out, index=False)

if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--genotypes", nargs="+")
    p.add_argument("--out")
    args = p.parse_args()
    main(args)