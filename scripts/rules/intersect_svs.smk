rule index_sv_vcf:
    input:
        vcf = "data/vcf/variants_freeze5_sv_insdel_sym.vcf.gz",
    output:
        "data/vcf/variants_freeze5_sv_insdel_sym.vcf.gz.tbi"
    shell:
        """
        tabix -p vcf {input.vcf}
        """

rule intersect_svs:
    input:
        vcf = "data/vcf/variants_freeze5_sv_insdel_sym.vcf.gz",
        vcf_idx = "data/vcf/variants_freeze5_sv_insdel_sym.vcf.gz.tbi",
        py_script = "scripts/find_sv_overlap.py",
        refseq = "data/mm39.refseq.tsv.gz"
    output:
        "csv/sv_gene_overlap.summary.csv"
    shell:
        """
        python {input.py_script} --refseq {input.refseq} \
                                 --vcf {input.vcf} \
                                 --out {output}
        """