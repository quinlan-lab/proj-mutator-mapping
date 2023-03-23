PROJDIR = "/Users/tomsasani/quinlanlab/proj-mutator-mapping"

BCFTOOLS = "/Users/tomsasani/bin/bin/bcftools"
TABIX = "/Users/tomsasani/bin/bin/tabix"

chroms_ = list(map(str, range(1, 20)))
chroms = ['chr' + c for c in chroms_]
 
crosses = ["bxd"]
samples = ["all", "conditioned"]
kmer_sizes = [1, 3]

rule all:
    input:
        #expand(PROJDIR + "/csv/{cross}/k{k}.genome.{samples}_samples.significant_markers.csv", cross=crosses, k=kmer_sizes, samples=samples,),
        PROJDIR + "/img/mgp_spectra_comparison.png"

rule download_cc_geno:
    input:
    output:
        temp(PROJDIR +"/data/genotypes/cc.{chrom}.geno")
    shell:
        """
        wget -O {output} https://raw.githubusercontent.com/rqtl/qtl2data/main/CC/cc_geno{wildcards.chrom}.csv
        """

rule combine_cc_geno:
    input:
        genotypes = expand(PROJDIR +"/data/genotypes/cc.{chrom}.geno", chrom=chroms_)
    output:
        out_geno = PROJDIR + "/data/genotypes/cc.geno"
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
        unformatted = PROJDIR + "/data/mutations/cc/cc.csv"
    output:
        formatted = PROJDIR + "/data/mutations/cc/annotated_filtered_singletons.csv"
    run:
        import pandas as pd 

        df = pd.read_csv(input.unformatted, dtype={"INDEL": int})
        df = df[df["INDEL"] != 1]

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
        orig = PROJDIR + "/data/genotypes/cc.markers.tmp"
    output:
        new = PROJDIR + "/data/genotypes/cc.markers"
    run:
        import pandas as pd 
        df = pd.read_csv(input.orig)
        df["Mb"] = df["position(b38)"] / 1_000_000
        df.to_csv(output.new, index=False)


rule download_singletons:
    input:
    output:
        temp(PROJDIR + "/data/mutations/{chrom}.singletons.csv")
    shell:
        """
        wget -O {output} https://raw.githubusercontent.com/tomsasani/bxd_mutator_manuscript/main/data/singleton_vars/{wildcards.chrom}_singleton_vars.exclude.csv
        """


rule combine_bxd_singletons:
    input:
        singletons = expand(PROJDIR + "/data/mutations/{chrom}.singletons.csv", chrom=chroms),
        py_script = PROJDIR + "/scripts/combine_bxd_singletons.py",
        metadata = PROJDIR + "/data/bam_names_to_metadata.xlsx",
        config = PROJDIR + "/data/json/{cross}.json",
    output:
        PROJDIR + "/data/mutations/{cross}/annotated_filtered_singletons.all_samples.csv"
    shell:
        """
        python {input.py_script} \
                        --metadata {input.metadata} \
                        --singletons {input.singletons} \
                        --config {input.config} \
                        --out {output}
        """

rule combine_bxd_singletons_conditional:
    input:
        singletons = expand(PROJDIR + "/data/mutations/{chrom}.singletons.csv", chrom=chroms),
        py_script = PROJDIR + "/scripts/combine_bxd_singletons.py",
        metadata = PROJDIR + "/data/bam_names_to_metadata.xlsx",
        config = PROJDIR + "/data/json/{cross}.json",
    output:
        PROJDIR + "/data/mutations/{cross}/annotated_filtered_singletons.conditioned_samples.csv"
    shell:
        """
        python {input.py_script} \
                        --metadata {input.metadata} \
                        --singletons {input.singletons} \
                        --config {input.config} \
                        --out {output} \
                        -condition_on_mutator
        """

rule run_ihd:
    input:
        singletons = PROJDIR + "/data/mutations/{cross}/annotated_filtered_singletons.{samples}_samples.csv",
        geno = PROJDIR + "/data/genotypes/{cross}.geno",
        config = PROJDIR + "/data/json/{cross}.json",
        py_script = PROJDIR + "/ihd/run_ihd_scan.py"
    output:
        PROJDIR + "/csv/{cross}.k{k}.genome.{samples}_samples.results.csv"
    shell:
        """
        python {input.py_script} --mutations {input.singletons} \
                                 --config {input.config} \
                                 --out {output} \
                                 -k {wildcards.k} \
                                 -permutations 1000 \
        """

rule plot_ihd:
    input:
        results = PROJDIR + "/csv/{cross}.k{k}.genome.{samples}_samples.results.csv",
        markers = PROJDIR + "/data/genotypes/{cross}.markers",
        py_script = PROJDIR + "/ihd/plot_ihd_results.py"
    output:
        PROJDIR + "/csv/{cross}/k{k}.genome.{samples}_samples.significant_markers.csv"
    shell:
        """
        python {input.py_script} --markers {input.markers} \
                                 --results {input.results} \
                                 --outpref /Users/tomsasani/quinlanlab/proj-mutator-mapping/csv/{wildcards.cross}/k{wildcards.k}.genome.{wildcards.samples}_samples \
        """

rule download_mgp_vcf:
    input:
        BCFTOOLS_PATH = BCFTOOLS,

    output:
        PROJDIR + "/data/vcf/mgp.regions.vcf.gz"
    shell:
        """
        {input.BCFTOOLS_PATH} view -Oz -o {output} -r 4:100000000-120000000,6:100000000-120000000 \
            ftp://ftp.ebi.ac.uk/pub/databases/mousegenomes/REL-1505-SNPs_Indels/mgp.v5.merged.snps_all.dbSNP142.vcf.gz
        """

rule tabix_mgp_vcf:
    input:
        vcf = PROJDIR + "/data/vcf/mgp.regions.vcf.gz",
        TABIX_PATH = TABIX,
    output:
        PROJDIR + "/data/vcf/mgp.regions.vcf.gz.tbi"
    shell:
        """
        {input.TABIX_PATH} -p vcf {input.vcf}
        """

rule compare_mgp_spectra:
    input:
        dumont_data = PROJDIR + "/data/SuppTables_concat.xlsx",
        mgp_vcf = PROJDIR + "/data/vcf/mgp.regions.vcf.gz",
        mgp_idx = PROJDIR + "/data/vcf/mgp.regions.vcf.gz.tbi",
        py_script = PROJDIR + "/scripts/make_tidy_df_mgp.py",
    output:
        PROJDIR + "/img/mgp_spectra_comparison.png"
    shell:
        """
        python {input.py_script} --dumont_xls {input.dumont_data} \
                                 --vcf {input.mgp_vcf} \
                                 --out {output}
        """