include: "rules/combine_singleton_data.smk"
include: "rules/run_ihd.smk"
include: "rules/process_mgp_data.smk"
include: "rules/process_bxd_data.smk"
include: "rules/calculate_callable_kmer.smk"

rule all:
    input:
        "figs/bxd/k1.genome.all_samples.png", # Figure 2a
        "figs/bxd/k1.genome.conditioned_samples.png", # Figure 2b
        "figs/bxd_spectra_1mer.png", # Figure 3a
        "figs/bxd_spectra_3mer.png", # Figure 3b
        "figs/bxd/k3.genome.all_samples.png", # Figure 2-supplement 1a
        "figs/bxd/k3.genome.conditioned_samples.png", # Figure 2-supplement 1b
        "figs/mgp_spectra_1mer.png", # Figure 3-supplement 1



