PROJDIR = "/scratch/ucgd/lustre-work/quinlan/u1006375/proj-mutator-mapping"
chroms = list(map(str, range(1, 23)))

rule all:
    input:
        expand(PROJDIR + "/csv/ceph/k{k}.genome.significant_markers.csv", k=[1,3]),

rule get_genos:
    input:
        py_script = PROJDIR + "/scripts/generate_ceph_geno.py"
    output:
        temp(PROJDIR + "/data/genotypes/per-chrom/{chrom}.ceph.geno"),
        temp(PROJDIR + "/data/genotypes/per-chrom/{chrom}.ceph.markers")
    shell:
        """
        python {input.py_script} --chrom {wildcards.chrom}
        """

rule combine_genos:
    input:
        genos = expand(PROJDIR + "/data/genotypes/per-chrom/{chrom}.ceph.geno", chrom=chroms),
    output:
        PROJDIR + "/data/genotypes/ceph.geno",
    shell:
        """
        for f in {input.genos}; do 
            head -n1 $f > g_header && tail -n+2 $f > g_footer && cat g_footer >> g_tmp_out; 
        done

        cat g_header g_tmp_out > {output}
        rm g_header
        rm g_footer
        rm g_tmp_out
        
        """

rule combine_markers:
    input:
        markers = expand(PROJDIR + "/data/genotypes/per-chrom/{chrom}.ceph.markers", chrom=chroms)
    output:
        PROJDIR + "/data/genotypes/ceph.markers"
    shell:
        """
        for f in {input.markers}; do 
            head -n1 $f > m_header && tail -n+2 $f > m_footer && cat m_footer >> m_tmp_out; 
        done

        cat m_header m_tmp_out > {output}
        rm m_header
        rm m_footer
        rm m_tmp_out
        
        """

rule run_manhattan:
    input:
        mutations = PROJDIR + "/data/mutations/ceph/formatted_dnms.csv",
        config = PROJDIR + "/data/json/ceph.json",
        py_script = PROJDIR + "/scripts/run_ihd_scan.py",
        genos = PROJDIR + "/data/genotypes/ceph.geno",
        markers = PROJDIR + "/data/genotypes/ceph.markers",
    output:
        PROJDIR + "/csv/ceph.k{k}.genome.results.csv"
    shell:
        """
        python {input.py_script} --mutations {input.mutations} \
                                 --config {input.config} \
                                 --out {output} \
                                 -k {wildcards.k} \
                                 -permutations 1000
        """

rule plot_manhattan:
    input:
        results = PROJDIR + "/csv/ceph.k{k}.genome.results.csv",
        markers = PROJDIR + "/data/genotypes/ceph.markers",
        py_script = PROJDIR + "/scripts/plot_ihd_results.py"
    output:
        PROJDIR + "/csv/ceph/k{k}.genome.significant_markers.csv"
    shell:
        """
        python {input.py_script} --markers {input.markers} \
                                 --results {input.results} \
                                 --outpref /scratch/ucgd/lustre-work/quinlan/u1006375/proj-mutator-mapping/csv/ceph/k{wildcards.k}.genome \
        """