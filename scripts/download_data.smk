PROJDIR = "/scratch/ucgd/lustre-work/quinlan/u1006375/proj-mutator-mapping"

chroms_ = list(map(str, range(1, 20)))
chroms = ['chr' + c for c in chroms_]

crosses = ["bxd"]

rule all:
    input:
        expand(PROJDIR + "/csv/{cross}/k{k}.genome.significant_markers.csv", k=[1,3], cross=crosses),
        #expand(PROJDIR + "/csv/per-chrom/{cross}/k{k}.{chrom}.significant_markers.csv", k=[1,3], chrom=chroms, cross=["bxd", "cc"]),


rule download_singletons:
    input:
    output:
        temp(PROJDIR + "/data/singletons/{chrom}.singletons.csv")
    shell:
        """
        wget -O {output} https://raw.githubusercontent.com/tomsasani/bxd_mutator_manuscript/main/data/singleton_vars/{wildcards.chrom}_singleton_vars.exclude.csv
        """

rule combine_bxd_singletons:
    input:
        singletons = expand(PROJDIR + "/data/singletons/{chrom}.singletons.csv", chrom=chroms),
        py_script = PROJDIR + "/scripts/combine_bxd_singletons.py",
        metadata = PROJDIR + "/data/bam_names_to_metadata.xlsx",
        config = PROJDIR + "/data/json/{cross}.json",
    output:
        PROJDIR + "/data/singletons/{cross}/annotated_filtered_singletons.csv"
    shell:
        """
        python {input.py_script} \
                        --metadata {input.metadata} \
                        --singletons {input.singletons} \
                        --config {input.config} \
                        --out {output}
        """

rule download_cc_geno:
    input:
    output:
        temp(PROJDIR + "/data/genotypes/cc/{chrom}.geno")
    shell:
        """
        wget -O {output} https://raw.githubusercontent.com/rqtl/qtl2data/main/CC/cc_geno{wildcards.chrom}.csv
        """

rule combine_cc_geno:
    input:
        genotypes = expand(PROJDIR + "/data/genotypes/cc/{chrom}.geno", chrom=chroms_),
        py_script = PROJDIR + "/scripts/combine_cc_genotypes.py"
    output:
        PROJDIR + "/data/genotypes/CC.geno"
    shell:
        """
        python {input.py_script} --genotypes {input.genotypes} --out {output}
        """

rule filter_cc_singletons:
    input:
        PROJDIR + "/data/singletons/cc/annotated_filtered_singletons.raw.csv"
    output:
        PROJDIR + "/data/singletons/cc/annotated_filtered_singletons.csv"
    shell:
        """
        grep -v "CC038\|CC062\|CC072" {input} > {output} 
        """

rule run_manhattan:
    input:
        singletons = PROJDIR + "/data/singletons/{cross}/annotated_filtered_singletons.csv",
        config = PROJDIR + "/data/json/{cross}.json",
        py_script = PROJDIR + "/scripts/compute_distance.py"
    output:
        PROJDIR + "/csv/{cross}.k{k}.genome.results.csv"
    shell:
        """
        python {input.py_script} --singletons {input.singletons} \
                                 --config {input.config} \
                                 --out {output} \
                                 -k {wildcards.k} \
        """

rule plot_manhattan:
    input:
        results = PROJDIR + "/csv/{cross}.k{k}.genome.results.csv",
        markers = PROJDIR + "/data/genotypes/{cross}.markers",
        py_script = PROJDIR + "/scripts/manhattan.py"
    output:
        PROJDIR + "/csv/{cross}/k{k}.genome.significant_markers.csv"
    shell:
        """
        python {input.py_script} --markers {input.markers} \
                                 --results {input.results} \
                                 --outpref /scratch/ucgd/lustre-work/quinlan/u1006375/proj-mutator-mapping/csv/{wildcards.cross}/k{wildcards.k}.genome \
        """

rule run_manhattan_chrom:
    input:
        singletons = PROJDIR + "/data/singletons/{cross}/annotated_filtered_singletons.csv",
        config = PROJDIR + "/data/json/{cross}.json",
        py_script = PROJDIR + "/scripts/compute_distance.py"
    output:
        PROJDIR + "/csv/per-chrom/{cross}/k{k}.{chrom}.results.csv"
    shell:
        """
        python {input.py_script} --singletons {input.singletons} \
                                 --config {input.config} \
                                 --out {output} \
                                 -k {wildcards.k} \
                                 -chrom {wildcards.chrom} \
                                 -permutations 1000
        """

rule plot_manhattan_chrom:
    input:
        results = PROJDIR + "/csv/per-chrom/{cross}/k{k}.{chrom}.results.csv",
        markers = PROJDIR + "/data/genotypes/{cross}.markers",
        py_script = PROJDIR + "/scripts/manhattan.py"
    output:
        PROJDIR + "/csv/per-chrom/{cross}/k{k}.{chrom}.significant_markers.csv"
    shell:
        """
        python {input.py_script} --markers {input.markers} \
                                 --results {input.results} \
                                 --outpref /scratch/ucgd/lustre-work/quinlan/u1006375/proj-mutator-mapping/csv/per-chrom/{wildcards.cross}/k{wildcards.k}.{wildcards.chrom} 
        """