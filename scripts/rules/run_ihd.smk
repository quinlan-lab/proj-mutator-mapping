rule run_ihd:
    input:
        singletons = "data/mutations/{cross}/annotated_filtered_singletons.{samples}_samples.csv",
        geno = "data/genotypes/{cross}.geno",
        config = "data/json/{cross}.json",
        py_script = "ihd/run_ihd_scan.py"
    output:
        "csv/{cross}.k{k}.genome.{samples}_samples.results.csv"
    shell:
        """
        python {input.py_script} --mutations {input.singletons} \
                                 --config {input.config} \
                                 --out {output} \
                                 -k {wildcards.k} \
                                 -distance_method cosine \
                                 -permutations 1000 \
                                 -progress \
                                 -threads 4
        """

rule plot_ihd:
    input:
        results = "csv/{cross}.k{k}.genome.{samples}_samples.results.csv",
        markers = "data/genotypes/{cross}.markers",
        py_script = "ihd/plot_ihd_results.py"
    output: "figs/{cross}.k{k}.genome.{samples}_samples.{ext}"
    shell:
        """
        python {input.py_script} --markers {input.markers} \
                                 --results {input.results} \
                                 --out {output} \
        """