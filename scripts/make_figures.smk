include: "rules/combine_singleton_data.smk"
include: "rules/run_ihd.smk"
include: "rules/process_mgp_data.smk"
include: "rules/process_bxd_data.smk"
include: "rules/calculate_callable_kmer.smk"
include: "rules/run_simulations.smk"
include: "rules/run_qtl_scans.smk"
include: "rules/run_sigprofiler.smk"

rule all:
    input:
        "figs/power_simulations.png", # Figure 1
        "figs/bxd.k1.genome.all_samples.eps", # Figure 2a
        "figs/bxd.k1.genome.conditioned_samples.eps", # Figure 2b

        "figs/bxd_spectra_1mer.C_A.eps", # Figure 3a
                "figs/bxd_spectra_1mer.C_A.png", # Figure 3a

        "figs/bxd_spectra_vs_age.png", # Figure 3b
        "figs/bxd_spectra_3mer_all.eps", # Figure 3c
        "figs/signature_activity.eps",
        "figs/bxd_spectra_1mer_all.eps", # Figure 3-supplement 1

        expand("figs/qtl_scans/{mutation_type}.png", mutation_type = ["C_T", "C_A", "C_G", "CpG_TpG", "A_T", "A_G", "A_C"]),

        "figs/mgp_spectra_1mer.eps", # Figure 3-supplement 2
        #"figs/power_comparison.TCC_TTC.png",
        "figs/power_comparison.png",



