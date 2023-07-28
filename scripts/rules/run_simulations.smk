import numpy as np
from skbio.stats.composition import clr 

# define parameter space for simulations
number_of_markers = [1_000]
number_of_haplotypes = [50]
number_of_mutations = [5, 20, 100]
number_of_permutations = [100]
mutation_types = ["C_T", "C_A"]
effect_sizes = list(range(100, 160, 10))
distance_methods = ["cosine"]#, "chisquare"]
linkages = [10, 25, 50, 100]
exp_afs = [50]

N_TRIALS = 100


rule run_ind_simulation:
    input:
        py_script = "ihd/run_ihd_power_simulations.py"
    output:
        temp("csv/power.{n_markers}.{n_haplotypes}.{n_mutations}.{n_permutations}.{effect_size}.{distance_method}.{mutation_type}.{linkage}.{exp_af}.csv"),
    shell:
        """
        python {input.py_script} --results {output} \
                                 -n_markers {wildcards.n_markers} \
                                 -n_haplotypes {wildcards.n_haplotypes} \
                                 -n_mutations {wildcards.n_mutations} \
                                 -n_permutations {wildcards.n_permutations} \
                                 -effect_size {wildcards.effect_size} \
                                 -distance_method {wildcards.distance_method} \
                                 -mutation_type {wildcards.mutation_type} \
                                 -tag_strength {wildcards.linkage} \
                                 -exp_af {wildcards.exp_af}
        """

rule combine_ind_simulations:
    input:
        results = expand("csv/power.{n_markers}.{n_haplotypes}.{n_mutations}.{n_permutations}.{effect_size}.{distance_method}.{mutation_type}.{linkage}.{exp_af}.csv", 
                n_markers = number_of_markers,
                n_haplotypes = number_of_haplotypes,
                n_mutations = number_of_mutations,
                n_permutations = number_of_permutations,
                effect_size = effect_sizes,
                distance_method = distance_methods,
                mutation_type = mutation_types,
                linkage = linkages,
                exp_af = exp_afs)
    output:
        combined_results = "csv/ihd_power.csv"
    run:
        import pandas as pd 
        dfs = []
        for fh in input.results:
            df = pd.read_csv(fh)
            dfs.append(df)
        dfs = pd.concat(dfs)
        dfs.to_csv(output.combined_results, index=False)


rule plot_simulations:
    input:
        py_script = "ihd/plot_power.py",
        results = "csv/ihd_power.csv"
    output: "figs/power_simulations.png"
    shell:
        """
        python {input.py_script} --results {input.results} \
                                 --out {output}
        """


rule generate_simulation_comp_data:
    input:
        py_script = "ihd/run_ihd_power_simulations.py"
    output:
        results = temp("data/Rqtl_sim/power_ind.{n_markers}.{n_haplotypes}.{n_mutations}.1.{effect_size}.{distance_method}.{mutation_type}.{linkage}.{exp_af}.csv"),
        genos = temp("data/Rqtl_sim/sim.{n_markers}.{n_haplotypes}.{n_mutations}.1.{effect_size}.{distance_method}.{mutation_type}.{linkage}.{exp_af}.geno"),
        spectra = temp("data/Rqtl_sim/spectra.{n_markers}.{n_haplotypes}.{n_mutations}.1.{effect_size}.{distance_method}.{mutation_type}.{linkage}.{exp_af}.csv")
    shell:
        """
        python {input.py_script} --results {output.results} \
                                 -n_markers {wildcards.n_markers} \
                                 -n_haplotypes {wildcards.n_haplotypes} \
                                 -n_mutations {wildcards.n_mutations} \
                                 -n_permutations 1 \
                                 -effect_size {wildcards.effect_size} \
                                 -distance_method {wildcards.distance_method} \
                                 -mutation_type {wildcards.mutation_type} \
                                 -tag_strength {wildcards.linkage} \
                                 -exp_af {wildcards.exp_af} \
                                 -raw_geno {output.genos} \
                                 -raw_spectra {output.spectra} \
        """


rule make_simulated_pmap:
    input:
    output: pmap = temp("data/Rqtl_sim/split/sim.{n_markers}.pmap")
    run:
        import numpy as np
        with open(output.pmap, "w") as outfh:
            header = ["marker", "chr", "pos"]
            print (",".join(header), file=outfh)
            positions = np.linspace(1, 200, num=int(wildcards.n_markers))
            for mi in range(int(wildcards.n_markers)):
                vals = [f"rs{mi}", "1", str(positions[mi])]
                print (",".join(vals), file=outfh)


rule make_simulated_gmap:
    input:
    output: gmap = temp("data/Rqtl_sim/split/sim.{n_markers}.gmap")
    run:
        import numpy as np
        n_markers = 100
        with open(output.gmap, "w") as outfh:
            header = ["marker", "chr", "pos"]
            print (",".join(header), file=outfh)
            positions = np.linspace(1, 100, num=int(wildcards.n_markers))
            for mi in range(int(wildcards.n_markers)):
                vals = [f"rs{mi}", "1", str(positions[mi])]
                print (",".join(vals), file=outfh)


rule reformat_simulated_spectra:
    input:
        spectra = "data/Rqtl_sim/spectra.{n_markers}.{n_haplotypes}.{n_mutations}.1.{effect_size}.{distance_method}.{mutation_type}.{linkage}.{exp_af}.csv"
    output:
        spectra = temp("data/Rqtl_sim/sim.{n_markers}.{n_haplotypes}.{n_mutations}.1.{effect_size}.{distance_method}.{mutation_type}.{linkage}.{exp_af}.pheno")
    run:
        import pandas as pd
        spectra_df = pd.read_csv(input.spectra)
        group_cols = ["sample", "trial", "focal_marker"]
        # convert spectra df to tidy format
        spectra_df_tidy = spectra_df.melt(id_vars=group_cols)
        # calculate total sums of each mutation type
        spectra_df_grouped = spectra_df_tidy.groupby(group_cols).agg({"value": sum}).reset_index().rename(columns={"value": "total"})
        spectra_df_combined = spectra_df_tidy.merge(spectra_df_grouped, on=group_cols)
        # convert counts to fractions
        spectra_df_combined["fraction"] = spectra_df_combined["value"] / spectra_df_combined["total"]
        # pivot df back to long-form
        spectra_df_long = spectra_df_combined.drop(columns=["total", "value"]).pivot(index=group_cols, columns="variable").reset_index()
        new_cols = ["sample", "trial", "focal_marker"]
        mut_cols =[c[1].replace(">", "_") for c in spectra_df_long.columns if c[0] == "fraction"]
        spectra_df_long.columns = new_cols + mut_cols

        spectra_df_long.to_csv(output.spectra, index=False)

rule split_spectra_by_trial:
    input:
        spectra = "data/Rqtl_sim/sim.{n_markers}.{n_haplotypes}.{n_mutations}.1.{effect_size}.{distance_method}.{mutation_type}.{linkage}.{exp_af}.pheno"
    output:
        spectra = temp("data/Rqtl_sim/split/sim.{n_markers}.{n_haplotypes}.{n_mutations}.1.{effect_size}.{distance_method}.{mutation_type}.{linkage}.{exp_af}.{trial}.pheno")
    run:
        import pandas as pd 

        df = pd.read_csv(input.spectra)
        df["is_in_trial"] = df["trial"].apply(lambda t: 1 if str(t) in list({wildcards.trial}) else 0)
        df = df.query('is_in_trial == 1')
        df = df.drop(columns=["trial", "is_in_trial"])
        df.to_csv(output.spectra, index=False)

rule split_geno_by_trial:
    input:
        geno = "data/Rqtl_sim/sim.{n_markers}.{n_haplotypes}.{n_mutations}.1.{effect_size}.{distance_method}.{mutation_type}.{linkage}.{exp_af}.geno"
    output:
        geno = temp("data/Rqtl_sim/split/sim.{n_markers}.{n_haplotypes}.{n_mutations}.1.{effect_size}.{distance_method}.{mutation_type}.{linkage}.{exp_af}.{trial}.geno")
    run:
        import pandas as pd 

        df = pd.read_csv(input.geno)
        df["is_in_trial"] = df["trial"].apply(lambda t: 1 if str(t) in list({wildcards.trial}) else 0)
        df = df.query('is_in_trial == 1')
        df = df.drop(columns=["trial", "is_in_trial"])
        df.to_csv(output.geno, index=False)



rule make_simulated_json:
    input:
        #spectra = "data/Rqtl_sim/sim.{n_markers}.{n_haplotypes}.{n_mutations}.1.{effect_size}.{distance_method}.{mutation_type}.{linkage}.{exp_af}.{trial}.pheno",
        geno = "data/Rqtl_sim/split/sim.{n_markers}.{n_haplotypes}.{n_mutations}.1.{effect_size}.{distance_method}.{mutation_type}.{linkage}.{exp_af}.{trial}.geno"
    output:
        out_json = temp("data/Rqtl_sim/split/sim.{n_markers}.{n_haplotypes}.{n_mutations}.1.{effect_size}.{distance_method}.{mutation_type}.{linkage}.{exp_af}.{trial}.json")
    run:
        import json 

        json_dict = {
            "description": "simulated BXD data",
            "crosstype": "risib",
            "sep": ",",
            "na.strings": ["-", "NA"],
            "comment.char": "#",
            "geno": f"sim.{wildcards.n_markers}.{wildcards.n_haplotypes}.{wildcards.n_mutations}.1.{wildcards.effect_size}.{wildcards.distance_method}.{wildcards.mutation_type}.{wildcards.linkage}.{wildcards.exp_af}.{wildcards.trial}.geno",
            "pmap": f"sim.{wildcards.n_markers}.pmap",
            "gmap": f"sim.{wildcards.n_markers}.gmap",
            "alleles": ["B", "D"],
            "genotypes": {"B": 1, "D": 2},
            "geno_transposed": True,
        }
        json_data = json.dumps(json_dict)
        json_file = open(output.out_json, "w")
        json_file.write(json_data)
        json_file.close()


rule map_simulated_qtl:
    input:
        genotypes = "data/Rqtl_sim/split/sim.{n_markers}.{n_haplotypes}.{n_mutations}.1.{effect_size}.{distance_method}.{mutation_type}.{linkage}.{exp_af}.{trial}.geno",
        spectra = "data/Rqtl_sim/split/sim.{n_markers}.{n_haplotypes}.{n_mutations}.1.{effect_size}.{distance_method}.{mutation_type}.{linkage}.{exp_af}.{trial}.pheno",
        pmap = "data/Rqtl_sim/split/sim.{n_markers}.pmap",
        gmap = "data/Rqtl_sim/split/sim.{n_markers}.gmap",
        json = "data/Rqtl_sim/split/sim.{n_markers}.{n_haplotypes}.{n_mutations}.1.{effect_size}.{distance_method}.{mutation_type}.{linkage}.{exp_af}.{trial}.json",
        Rscript = "scripts/map_qtl.R",
    output:
        power = temp("csv/simulation_results_rqtl.{n_markers}.{n_haplotypes}.{n_mutations}.1.{effect_size}.{distance_method}.{mutation_type}.{linkage}.{exp_af}.{trial}.csv")
    shell: 
        """
        Rscript {input.Rscript} -j {input.json} -p {input.spectra} -o {output.power} -m {wildcards.mutation_type} -t {wildcards.trial}
        """

rule combine_rqtl_simulation_results:
    input: 
        power = expand("csv/simulation_results_rqtl.{n_markers}.{n_haplotypes}.{n_mutations}.1.{effect_size}.{distance_method}.{mutation_type}.{linkage}.{exp_af}.{trial}.csv", 
                n_markers = number_of_markers,
                n_haplotypes = number_of_haplotypes,
                n_mutations = number_of_mutations,
                effect_size = effect_sizes,
                distance_method = distance_methods,
                mutation_type = mutation_types,
                linkage = linkages,
                exp_af = exp_afs,
                trial = list(range(N_TRIALS)))
    output:
        results = "csv/rqtl_power.csv"
    run:
        dfs = []
        for fh in input.power:
            (n_markers, n_haplotypes, n_mutations, n_permutations, effect_size, distance_method, mutation_type, linkage, exp_af, trial) = fh.split('/')[-1].split('.')[1:-1]
            df = pd.read_csv(fh)
            df["n_markers"] = n_markers 
            df["n_haplotypes"] = n_haplotypes
            df["n_mutations"] = n_mutations 
            df["n_permutations"] = n_permutations 
            df["effect_size"] = effect_size 
            df["distance_method"] = distance_method 
            df["mutation_type"] = mutation_type
            df["tag_strength"] = linkage
            df["exp_af"] = exp_af
            df["trial"] = trial
            dfs.append(df)
        dfs = pd.concat(dfs)
        dfs.to_csv(output.results, index=False)


rule plot_power_comparison:
    input:
        ihd_power = "csv/ihd_power.csv",
        qtl_power = "csv/rqtl_power.csv",
        py_script = "scripts/compare_ihd_qtl_power.py"
    output: "figs/power_comparison.png"
    shell:
        """
        python {input.py_script} --qtl_power {input.qtl_power} --ihd_power {input.ihd_power} --out {output}
        """


rule plot_cosine_vs_chisquare:
    input:
        ihd_power = "csv/ihd_power.csv",
        py_script = "scripts/compare_cosine_chisquare_power.py"
    output: "figs/power_comparison_distance_methods.png"
    shell:
        """
        python {input.py_script} --ihd_power {input.ihd_power} --out {output}
        """