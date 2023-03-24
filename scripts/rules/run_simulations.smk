rule run_simulations:
    input:
        py_script = "ihd/run_ihd_power_simulations.py"
    output:
        "csv/simulated_power.csv"
    shell:
        """
        python {input.py_script} --out {output}
        """

rule plot_simulations:
    input:
        py_script = "ihd/plot_power.py",
        results = "csv/simulated_power.csv"
    output:
        "figs/power_simulations.png"
    shell:
        """
        python {input.py_script} --results {input.results} \
                                 --out {output}
        """