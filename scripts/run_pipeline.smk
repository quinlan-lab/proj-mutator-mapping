PROJDIR = "/Users/tomsasani/quinlanlab/proj-mutator-mapping"

chroms_ = list(map(str, range(1, 20)))
chroms = ['chr' + c for c in chroms_]
 
rule all:
    input:
        expand(PROJDIR + "/csv/bxd/k{k}.genome.significant_markers.csv", k=[1]),
        PROJDIR + "/data/mutations/bxd/annotated_filtered_shared.csv"

rule download_singletons:
    input:
    output:
        temp(PROJDIR + "/data/mutations/{chrom}.singletons.csv")
    shell:
        """
        wget -O {output} https://raw.githubusercontent.com/tomsasani/bxd_mutator_manuscript/main/data/singleton_vars/{wildcards.chrom}_singleton_vars.exclude.csv
        """

rule download_shared:
    input:
    output:
        temp(PROJDIR + "/data/mutations/{chrom}.shared.csv")
    shell:
        """
        wget -O {output} https://raw.githubusercontent.com/tomsasani/bxd_mutator_manuscript/main/data/shared_vars/{wildcards.chrom}_shared_vars.exclude.csv
        """

rule combine_bxd_singletons:
    input:
        singletons = expand(PROJDIR + "/data/mutations/{chrom}.singletons.csv", chrom=chroms),
        py_script = PROJDIR + "/scripts/combine_bxd_singletons.py",
        metadata = PROJDIR + "/data/bam_names_to_metadata.xlsx",
        config = PROJDIR + "/data/json/{cross}.json",
    output:
        PROJDIR + "/data/mutations/{cross}/annotated_filtered_singletons.csv"
    shell:
        """
        python {input.py_script} \
                        --metadata {input.metadata} \
                        --singletons {input.singletons} \
                        --config {input.config} \
                        --out {output}
        """

rule combine_bxd_shared:
    input:
        singletons = expand(PROJDIR + "/data/mutations/{chrom}.shared.csv", chrom=chroms),
        py_script = PROJDIR + "/scripts/combine_bxd_singletons.py",
        metadata = PROJDIR + "/data/bam_names_to_metadata.xlsx",
        config = PROJDIR + "/data/json/{cross}.json",
    output:
        PROJDIR + "/data/mutations/{cross}/annotated_filtered_shared.csv"
    shell:
        """
        python {input.py_script} \
                        --metadata {input.metadata} \
                        --singletons {input.singletons} \
                        --config {input.config} \
                        --out {output}
        """

rule run_manhattan:
    input:
        singletons = PROJDIR + "/data/mutations/{cross}/annotated_filtered_singletons.csv",
        config = PROJDIR + "/data/json/{cross}.json",
        py_script = PROJDIR + "/ihd/run_ihd_scan.py"
    output:
        PROJDIR + "/csv/{cross}.k{k}.genome.results.csv"
    shell:
        """
        python {input.py_script} --mutations {input.singletons} \
                                 --config {input.config} \
                                 --out {output} \
                                 -k {wildcards.k} \
                                 -permutations 100 
        """

rule plot_manhattan:
    input:
        results = PROJDIR + "/csv/{cross}.k{k}.genome.results.csv",
        markers = PROJDIR + "/data/genotypes/{cross}.markers",
        py_script = PROJDIR + "/ihd/plot_ihd_results.py"
    output:
        PROJDIR + "/csv/{cross}/k{k}.genome.significant_markers.csv"
    shell:
        """
        python {input.py_script} --markers {input.markers} \
                                 --results {input.results} \
                                 --outpref /Users/tomsasani/quinlanlab/proj-mutator-mapping/csv/{wildcards.cross}/k{wildcards.k}.genome \
        """