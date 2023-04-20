rule make_tidy_bxd_spectra_df:
    input:
        kmer_content = "data/combined.callable_kmer_content.csv",
        genos = "data/genotypes/bxd.geno",
        mutations = "data/mutations/bxd/annotated_filtered_singletons.condition_on_N.csv",
        py_script = "scripts/reformat_bxd_spectra.py"
    output: "csv/bxd_spectra_{k}mer.csv"
    shell:
        """
        python {input.py_script} --kmer_content {input.kmer_content} \
                                 --genos {input.genos} \
                                 --mutations {input.mutations} \
                                 --out {output} \
                                 -k {wildcards.k}
        """


rule plot_bxd_spectra_ca:
    input:
        spectra = "csv/bxd_spectra_{k}mer.csv",
        py_script = "scripts/plot_bxd_spectra.py"
    output: "figs/bxd_spectra_{k}mer.{mutation_type}.{ext}"
    shell:
        """
        python {input.py_script} --spectra {input.spectra} \
                                 --out {output} \
                                 -k {wildcards.k} \
                                 -mutation_type {wildcards.mutation_type} \
                                 -phenotype Rate
        """

rule plot_bxd_spectra_all:
    input:
        spectra = "csv/bxd_spectra_{k}mer.csv",
        py_script = "scripts/plot_bxd_spectra.py"
    output: "figs/bxd_spectra_{k}mer_all.{ext}"
    shell:
        """
        python {input.py_script} --spectra {input.spectra} \
                                 --out {output} \
                                 -k {wildcards.k} \
                                 -phenotype Fraction
        """

rule plot_bxd_spectra_vs_age:
    input:
        spectra = "csv/bxd_spectra_1mer.csv",
        py_script = "scripts/plot_bxd_spectra_vs_age.py"
    output: "figs/bxd_spectra_vs_age.{ext}"
    shell:
        """
        python {input.py_script} --spectra {input.spectra} \
                                 --out {output} \
                                 --mutation_type C_A \
                                 -phenotype Count
        """