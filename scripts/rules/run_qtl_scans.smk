rule run_qtl_scan:
    input:
        json = "data/Rqtl_data/bxd.rqtl2.json",
        spectra = "csv/bxd_spectra_1mer.csv",
        Rscript = "scripts/map_bxd_qtl.R"
    output: "figs/qtl_scans/{mutation_type}.png"
    shell:
        """
        Rscript {input.Rscript} -j {input.json} -m {wildcards.mutation_type} -p {input.spectra} -o {output}
        """