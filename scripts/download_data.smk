PROJDIR = "/scratch/ucgd/lustre-work/quinlan/u1006375/proj-mutator-mapping"

chroms_ = list(map(str, range(1, 20)))
chroms = ['chr' + c for c in chroms_]

rule all:
    input:
        expand(PROJDIR + "/csv/bxd.k{k}.significant_markers.csv", k=[1,3])

rule download_singletons:
    input:
    output:
        temp(PROJDIR + "/data/singletons/{chrom}.singletons.csv")
    shell:
        """
        wget -O {output} https://raw.githubusercontent.com/tomsasani/bxd_mutator_manuscript/main/data/singleton_vars/{wildcards.chrom}_singleton_vars.exclude.csv
        """

rule combine_singletons:
    input:
        singletons = expand(PROJDIR + "/data/singletons/{chrom}.singletons.csv", chrom=chroms),
        py_script = PROJDIR + "/scripts/combine_bxd_singletons.py",
        metadata = PROJDIR + "/data/bam_names_to_metadata.xlsx",
    output:
        PROJDIR + "/data/singletons/bxd/annotated_filtered_singletons.csv"
    shell:
        """
        python {input.py_script} \
                        --metadata {input.metadata} \
                        --singletons {input.singletons} \
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

rule run_manhattan:
    input:
        singletons = PROJDIR + "/data/singletons/bxd/annotated_filtered_singletons.csv",
        config = PROJDIR + "/data/json/bxd.json",
        py_script = PROJDIR + "/scripts/compute_distance.py"
    output:
        PROJDIR + "/csv/bxd.k{k}.results.csv"
    shell:
        """
        python {input.py_script} --singletons {input.singletons} \
                                 --config {input.config} \
                                 --out {output} \
                                 -k {wildcards.k} \
        """

rule plot_manhattan:
    input:
        results = PROJDIR + "/csv/bxd.k{k}.results.csv",
        markers = PROJDIR + "/data/genotypes/BXD.markers",
        py_script = PROJDIR + "/scripts/manhattan.py"
    output:
        PROJDIR + "/csv/bxd.k{k}.significant_markers.csv"
    shell:
        """
        python {input.py_script} --markers {input.markers} \
                                 --results {input.results} \
                                 --outpref /scratch/ucgd/lustre-work/quinlan/u1006375/proj-mutator-mapping/csv/bxd.k{wildcards.k} \
        """

