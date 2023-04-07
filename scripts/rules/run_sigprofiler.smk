rule reformat_spectra_for_spe:
    input:
        mutations = "data/mutations/bxd/annotated_filtered_singletons.condition_on_N.csv",
        py_script = "scripts/reformat_spectra_for_sigprofiler.py"
    output:
        "data/mutations/bxd/annotated_filtered_singletons.condition_on_N.sig_profiler.txt"
    shell: 
        """
        python {input.py_script} --mutations {input.mutations} --out {output}
        """

rule run_sigprofiler:
    input:
        mutations = "data/mutations/bxd/annotated_filtered_singletons.condition_on_N.sig_profiler.txt",
        py_script = "scripts/run_sigprofiler.py"
    output:
        "data/sig_profiler_results/SBS96/Suggested_Solution/COSMIC_SBS96_Decomposed_Solution/Activities/COSMIC_SBS96_Activities.txt"
    shell:
        """
        python {input.py_script} --spectra {input.mutations} --outdir data/sig_profiler_results/
        """

rule plot_sbs_activities:
    input:
        activities = "data/sig_profiler_results/SBS96/Suggested_Solution/COSMIC_SBS96_Decomposed_Solution/Activities/COSMIC_SBS96_Activities.txt",
        mutations = "csv/bxd_spectra_1mer.csv",
        py_script = "scripts/plot_signature_activity.py"
    output: "figs/signature_activity.{ext}"
    shell:
        """
        python {input.py_script} --spectra {input.mutations} \
                                 --activities {input.activities} \
                                 --out {output}
        """
    