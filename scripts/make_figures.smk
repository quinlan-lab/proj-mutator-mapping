include: "rules/combine_singleton_data.smk"
include: "rules/run_ihd.smk"
include: "rules/process_mgp_data.smk"
include: "rules/process_bxd_data.smk"
include: "rules/calculate_callable_kmer.smk"
include: "rules/run_simulations.smk"
include: "rules/run_qtl_scans.smk"
include: "rules/run_sigprofiler.smk"
include: "rules/intersect_svs.smk"

mutations = ["C_T", "C_A", "C_G", "CpG_TpG", "A_T", "A_G", "A_C"]

###
# Rule below generates all figures/CSV output using the rules 
# that are imported above!
###

rule all:
    input:
        # FIGURE 1
        "figs/power_simulations.png", # supplement 1
        "figs/power_comparison.png", # supplement 2 
        "figs/power_comparison_distance_methods.png", # supplement 3
        
        # FIGURE 2
        "figs/bxd.k1.genome.condition_on_N.eps", # panel A
        "figs/bxd.k1.genome.condition_on_D.eps", # panel B

        "figs/bxd.k1.genome.condition_on_B.eps", # supplement 1
        expand("figs/qtl_scans/{mutation_type}.png", 
                    mutation_type = mutations),  # supplement 2

        # FIGURE 3
        "figs/bxd_spectra_1mer.C_A.png", # panel A
        "figs/bxd_spectra_vs_age.png", # panel B
        "figs/signature_activity.eps", # panel C

        "figs/bxd_spectra_1mer_all.eps", # supplement 1
        "figs/mgp_spectra_1mer.eps",     # supplement 2

        # ANALYSIS OUTPUTS
        "csv/sv_gene_overlap.summary.csv" # SV overlaps with genes




