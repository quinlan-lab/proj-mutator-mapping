# Find fixed SVs between DBA/2J and C57BL/6J that might
# explain interesting mutator signal

PROJDIR = "/scratch/ucgd/lustre-work/quinlan/u1006375/proj-mutator-mapping"

rule all:
    input: 
        PROJDIR + "/data/vcf/combined/merged.regenotyped.vcf"

rule download_c57:
    input:
    output:
        PROJDIR + "/data/bam/C57BL_6NJ.bam"
    threads: 8
    shell:
        """
        samtools view -b -o {output} -@ 8 ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR988/ERR9880493/C57BL_6NJ.bam
        """

rule download_dba:
    input:
    output:
        PROJDIR + "/data/bam/DBA_2J.bam"
    threads: 8
    shell:
        """
        samtools view -b -o {output} -@ 8 ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR988/ERR9880647/DBA_2J.bam
        """

rule index_bam:
    input:
        PROJDIR + "/data/bam/{sample}.bam"
    output:
        PROJDIR + "/data/bam/{sample}.bam.bai"
    threads: 8 
    shell:
        """
        samtools index -@ 8 {input}
        """

rule sample_call:
    input:
        bam = PROJDIR + "/data/bam/{sample}.bam",
        bai = PROJDIR + "/data/bam/{sample}.bam.bai",
        ref = "/scratch/ucgd/lustre-work/quinlan/u1006375/proj-mouse-mutators/data/ref/mm39.fa",
    output:
        PROJDIR + "/data/vcf/{sample}.vcf"
    threads: 8
    resources:
        mem_mb = 32000,
        #disk_mb = 64000,
        runtime = 1200
    shell:
        """
        dysgu run -p 8 \
                  -v2 \
                  -x \
                  --min-support 5 \
                  --remap False \
                  {input.ref} \
                  /scratch/ucgd/lustre-work/quinlan/u1006375/dysgu_tmp/{wildcards.sample}_tmp/ \
                  {input.bam} > {output}
        """

rule sample_merge:
    input:
       vcfs = expand(PROJDIR + "/data/vcf/{sample}.vcf", sample = ["C57BL_6NJ", "DBA_2J"]),
    output:
        PROJDIR + "/data/vcf/combined/merged.vcf"
    shell:
        """
        dysgu merge {input.vcfs} > {output}
        """

rule re_genotype:
    input:
        bam = PROJDIR + "/data/bam/{sample}.bam",
        bai = PROJDIR + "/data/bam/{sample}.bam.bai",
        ref = "/scratch/ucgd/lustre-work/quinlan/u1006375/proj-mouse-mutators/data/ref/mm39.fa",
        sites = PROJDIR + "/data/vcf/combined/merged.vcf"
    output:
        PROJDIR + "/data/vcf/re-genotyped/{sample}.vcf"
    threads: 8
    resources:
        mem_mb = 32000,
        #disk_mb = 64000,
        runtime = 1200
    shell:
        """
        dysgu run -p 8 \
                  --sites {input.sites} \
                  -v2 \
                  -x \
                  --min-support 5 \
                  --remap False \
                  {input.ref} \
                  /scratch/ucgd/lustre-work/quinlan/u1006375/dysgu_tmp/{wildcards.sample}_regeno_tmp/ \
                  {input.bam} > {output}
        """

rule re_merge:
    input:
        vcfs = expand(PROJDIR + "/data/vcf/re-genotyped/{sample}.vcf", sample = ["C57BL_6NJ", "DBA_2J"]),
    output:
        PROJDIR + "/data/vcf/combined/merged.regenotyped.vcf"
    shell:
        """
        dysgu merge {input.vcfs} > {output}
        """

