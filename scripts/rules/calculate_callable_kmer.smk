import pandas as pd

meta = pd.read_excel("data/bam_names_to_metadata.xlsx").dropna()
samples = meta["bam_name"].unique()

rule download_ref:
    input:
    output:
        temp("data/ref/mm10.fa.gz")
    shell:
        """
        wget -O {output} \
        http://hgdownload.soe.ucsc.edu/goldenPath/mm10/bigZips/mm10.fa.gz
        """

rule unzip_ref:
    input:
        "data/ref/mm10.fa.gz"
    output:
        temp("data/ref/mm10.fa")
    shell:
        """
        gzip -d {input}
        """


rule download_exclude:
    input:
    output:
        "data/exclude/mm10.seg_dups.simple_repeats.merged.bed.gz"
    shell:
        """
        wget -O {output} \
        https://github.com/tomsasani/bxd_mutator_manuscript/blob/main/data/mm10.seg_dups.simple_repeats.merged.bed.gz?raw=true
        """

rule download_coverage:
    input:
    output:
        temp("data/coverage/{sample}.bam.thresholds.bed")
    shell:
        """
        wget -O {output} https://raw.githubusercontent.com/tomsasani/bxd_mutator_manuscript/main/data/coverage_files/{wildcards.sample}.bam.thresholds.bed
        """

rule calculate_kmer_content:
    input:
        coverage = "data/coverage/{sample}.bam.thresholds.bed",
        ref = "data/ref/mm10.fa",
        exclude = "data/exclude/mm10.seg_dups.simple_repeats.merged.bed.gz",
        script = "scripts/calculate_callable_kmer_content.py"
    output:
        temp("data/coverage/per-sample/{sample}.callable_kmer_content.csv")
    shell:
        """
        python {input.script} --coverage {input.coverage} \
                              --ref {input.ref} \
                              --exclude {input.exclude} \
                              --out {output}
        """

rule combine_kmer_content:
    input:
        fhs = expand("data/coverage/per-sample/{sample}.callable_kmer_content.csv", sample=samples),
        metadata = "data/bam_names_to_metadata.xlsx",
    output:
        kmer_content = "data/combined.callable_kmer_content.csv"
    run:
        res = []
        for fh in input.fhs:
            df = pd.read_csv(fh)
            res.append(df)
        res_df = pd.concat(res)

        metadata = pd.read_excel(input.metadata).dropna()

        res_df_merged = res_df.merge(
            metadata,
            left_on="sample",
            right_on="bam_name",
        )

        res_df_merged.to_csv(output.kmer_content, index=False)

    