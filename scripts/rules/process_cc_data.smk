rule download_cc_geno:
    input:
    output:
        temp("data/genotypes/cc.{chrom}.geno")
    shell:
        """
        wget -O {output} https://raw.githubusercontent.com/rqtl/qtl2data/main/CC/cc_geno{wildcards.chrom}.csv
        """

rule combine_cc_geno:
    input:
        genotypes = expand("data/genotypes/cc.{chrom}.geno", chrom=chroms_)
    output:
        out_geno = "data/genotypes/cc.geno"
    run:
        import pandas as pd 

        dfs = []
        for fh in input.genotypes:
            df = pd.read_csv(fh, skiprows=3)
            dfs.append(df)
        dfs = pd.concat(dfs)
        sample_names = dfs.columns[1:]
        new_names = ["marker"] + [c.split("/")[0] for c in sample_names]
        dfs.columns = new_names
        dfs.to_csv(output.out_geno, index=False)

rule format_cc_mutations:
    input:
        unformatted = "data/mutations/cc/Private_Variants.csv"
    output:
        formatted = "data/mutations/cc/annotated_filtered_singletons.condition_on_N.csv"
    run:
        import pandas as pd 

        df = pd.read_csv(input.unformatted, dtype={"INDEL": int})
        df = df[df["INDEL"] != 1]

        print (df)

        revcomp = {"T": "A", "A": "T", "C": "G", "G": "C"}
        base_nucs = ["A", "C"]

        def convert_to_kmer(row):
            ref, alt = row["REF"], row["ALT"]
            if ref not in base_nucs:
                ref, alt = revcomp[ref], revcomp[alt]
            return f"N{ref}N>N{alt}N"

        df["kmer"] = df.apply(lambda row: convert_to_kmer(row), axis=1)
        df.rename(columns={"ANIMAL_NAME": "sample"}, inplace=True)

        df = df[~df["sample"].isin(["CC062", "CC038", "CC072"])]

        df["count"] = 1
        df[["sample", "kmer", "count"]].to_csv(output.formatted, index=False)


rule fix_cc_markers:
    input:
        orig = "data/genotypes/cc.markers.tmp"
    output:
        new = "data/genotypes/cc.markers"
    run:
        import pandas as pd 
        df = pd.read_csv(input.orig)
        df["Mb"] = df["position(b38)"] / 1_000_000
        df.to_csv(output.new, index=False)