rule run_ihd:
    input:
        singletons = "data/mutations/{cross}/annotated_filtered_singletons.condition_on_{condition}.csv",
        geno = "data/genotypes/{cross}.geno",
        config = "data/json/{cross}.json",
        py_script = "ihd/run_ihd_scan.py"
    output: "csv/{cross}.k{k}.genome.condition_on_{condition}.results.csv"
    shell:
        """
        python {input.py_script} --mutations {input.singletons} \
                                 --config {input.config} \
                                 --out {output} \
                                 -k {wildcards.k} \
                                 -distance_method cosine \
                                 -permutations 10000 \
                                 -progress \
                                 -threads 4 \
                                 -stratify_column true_epoch 
        """

rule plot_ihd:
    input:
        results = "csv/{cross}.k{k}.genome.condition_on_{condition}.results.csv",
        markers = "data/genotypes/{cross}.markers",
        py_script = "ihd/plot_ihd_results.py"
    output: "figs/{cross}.k{k}.genome.condition_on_{condition}.{ext}"
    shell:
        """
        python {input.py_script} --markers {input.markers} \
                                 --results {input.results} \
                                 --out {output} \
                                 -scale 1
        """