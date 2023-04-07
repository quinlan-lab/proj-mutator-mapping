chroms_ = list(map(str, range(1, 20)))
chroms = ['chr' + c for c in chroms_]

rule download_singletons:
    input:
    output:
        temp("data/mutations/{chrom}.singletons.csv")
    shell:
        """
        wget -O {output} https://raw.githubusercontent.com/tomsasani/bxd_mutator_manuscript/main/data/singleton_vars/{wildcards.chrom}_singleton_vars.exclude.csv
        """


rule combine_bxd_singletons:
    input:
        singletons = expand("data/mutations/{chrom}.singletons.csv", chrom=chroms),
        py_script = "scripts/combine_bxd_singletons.py",
        metadata = "data/bam_names_to_metadata.xlsx",
        config = "data/json/{cross}.json",
    output: "data/mutations/{cross}/annotated_filtered_singletons.condition_on_{condition}.csv"
    shell:
        """
        python {input.py_script} \
                        --metadata {input.metadata} \
                        --singletons {input.singletons} \
                        --config {input.config} \
                        --out {output} \
                        -condition_on_mutator {wildcards.condition}
        """