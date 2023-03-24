rule make_tidy_bxd_spectra_df:
    input:
        kmer_content = "data/combined.callable_kmer_content.csv",
        genos = "data/genotypes/bxd.geno",
        mutations = "data/mutations/bxd/annotated_filtered_singletons.all_samples.csv",
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


rule plot_bxd_spectra:
    input:
        spectra = "csv/bxd_spectra_{k}mer.csv",
        py_script = "scripts/plot_bxd_spectra.py"
    output: "figs/bxd_spectra_{k}mer.png"
    shell:
        """
        python {input.py_script} --spectra {input.spectra} \
                                 --out {output} \
                                 -k {wildcards.k}
        """
