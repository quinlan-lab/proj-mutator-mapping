rule download_mgp_vcf:
    input:
    output:
        "data/vcf/mgp.regions.vcf.gz"
    shell:
        """
        bcftools view -Oz -o {output} -r 4:100000000-120000000,6:100000000-120000000 \
            ftp://ftp.ebi.ac.uk/pub/databases/mousegenomes/REL-1505-SNPs_Indels/mgp.v5.merged.snps_all.dbSNP142.vcf.gz
        """

rule tabix_mgp_vcf:
    input:
        vcf = "data/vcf/mgp.regions.vcf.gz",
    output:
        "data/vcf/mgp.regions.vcf.gz.tbi"
    shell:
        """
        tabix -p vcf {input.vcf}
        """

rule reformat_mgp_spectra:
    input:
        dumont_data = "data/SuppTables_concat.xlsx",
        mgp_vcf = "data/vcf/mgp.regions.vcf.gz",
        mgp_idx = "data/vcf/mgp.regions.vcf.gz.tbi",
        py_script = "scripts/reformat_mgp_spectra.py",
    output:
        "csv/mgp_spectra.csv"
    shell:
        """
        python {input.py_script} --dumont_xls {input.dumont_data} \
                                 --vcf {input.mgp_vcf} \
                                 --out {output}
        """

rule plot_mgp_spectra:
    input:
        spectra = "csv/mgp_spectra.csv",
        py_script = "scripts/plot_mgp_spectra.py"
    output: "figs/mgp_spectra_1mer.{ext}"
    shell:
        """
        python {input.py_script} --spectra {input.spectra} \
                                 --out {output}
        """