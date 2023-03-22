# Mapping mutator alleles with inter-haplotype distance

[![docs](https://img.shields.io/badge/docs-latest-blue.svg)](https://quinlan-lab.github.io/proj-mutator-mapping/reference/) 
![pytest](https://github.com/quinlan-lab/proj-mutator-mapping/actions/workflows/tests.yaml/badge.svg)

## Summary

Identify alleles that affect the mutation spectrum in bi-parental recombinant inbred crosses. 

### Overview of method

![](img/fig-distance-method.png)

> **Overview of inter-haplotype distance method.**
> Figure 4: Overview of inter-haplotype distance method for discovering mutator alleles. **a)** A population of four haplotypes has been genotyped at three informative markers ($g_1$ through $g_3$); each haplotype also harbors private *de novo* germline mutations, which are denoted as triangles. *De novo* mutations are further classified into three possible mutation types, indicated by the color of each triangle. At each informative marker $g_n$, we calculate the total number of each mutation type observed on haplotypes that carry either parental allele (i.e., the aggregate mutation spectrum). We then calculate a $\chi^2$ statistic between the two aggregate mutation spectra. We repeat this process for every informative marker $g_n$. **b)** To assess the significance of any $\chi^2$ statistic peaks in a), we perform a permutation test by shuffling the haplotype labels associated with the *de novo* mutation data. In each of $N$ permutations, we record the maximum $\chi^2$ statistic encountered at any locus in the distance scan. Finally, we calculate the $1 - p$ percentile of the distribution of those maximum statistics to obtain a genome-wide $\chi^2$ statistic threshold at the specified value of $p$.
## Requirements

### Python dependencies

These methods were written in Python 3.9.6, and the full list of dependencies is provided in `requirements.txt`.

Dependencies can be installed with `pip install -r requirements.txt`.

> I recommend using [`pyenv`](https://github.com/pyenv/pyenv) and [`pyenv-virtualenv`](https://github.com/pyenv/pyenv-virtualenv) to manage Python environments.

### Description of input files 

Before running an inter-haplotype distance (IHD) scan, you'll need to prepare a
small number of input files.

1. ***De novo* germline mutation data**

    Mutation data should be in a CSV file with **three required columns** as follows:

    | sample | kmer | count |
    | - | - | - |
    | sample_A | CCT>CAT | 1 |
    | sample_A | TGG>TCG | 1 |
    | sample_B | GCA>GAA | 1 |


    **Notes:**

    > The `kmer` column *must* contain mutation types in the 3-mer format shown above -- the format will be validated at runtime.

    > The CSV file can contain either a) one row for every individual mutation observed in each sample, in which case the `count` column should always be set to 1 or b) the aggregate count of every mutation type observed in the sample, in which case the `count` column will reflect the total number of each mutation type observed in the sample.

    > The dataframe can have any number of additional columns, but only the three defined above will be used.

2. **Marker genotypes**

    Genotype data should be formatted in a similar fashion as in [R/qtl2](https://kbroman.org/qtl2/). Genotypes should be in a CSV file with N rows, where N is the number of genotyped markers. There should be a single column denoting the marker name, and as many columns as there are samples. See below:

    | marker | sample_A | sample_B | sample_C |
    | - | - | - | - |
    | rs001 | B | B | D |
    | rs002 | H | D | B |
    | rs003 | D | D | D |


3. **Marker information (optional)**

    If you wish to generate Manhattan-esque plots that summarize the results
    of an IHD scan, you'll need to provide a final CSV that links marker IDs with
    either physical or genetic map positions (or both). This file should contain a column called `marker`, a column called `chromosome`, and a column specifying one or both of `cM` or `Mb`.

    | marker | chromosome | cM | Mb |
    | - | - | - | - |
    | rs001 | 1 | 1 | 4.5230 |
    | rs002 | 1 | 2.6 | 5.1994 |
    | rs003 | 1 | 2.8 | 5.4872 |


4. **Configuration file**

    The final required file specifies the *absolute* paths to the marker genotypes and associated metadata, and defines the genotypes that are present in the marker genotype file. This configuration file should be JSON-formatted and look something like this:

    ```
    {
        "genotypes": {
            "B": 0,
            "D": 2,
            "H": 1
        },
        "geno": "path/to/geno/csv",
        "markers": "path/to/marker/metadata/csv"
    }
    ```

    **Notes:**

    > The `genotypes` dictionary should map the observed genotypes in file #2 to integer values that will be used during the IHD scan.

    > The two parental alleles *must* be mapped to values of 0 and 2, respectively. Heterozygous and unknown genotypes *must* be mapped to values of 1.

## Usage

### Running the full pipeline on BXD data

Using *de novo* germline mutation data from the BXD recombinant inbred mouse lines (originally generated in [Sasani et al. [2022]](https://www.nature.com/articles/s41586-022-04701-5)), the IHD scan and plotting scripts can be run in a single command using `snakemake` as follows:

```
snakemake -j1 -s scripts/run_pipeline.smk
```

If desired, the  `-j` parameter can be used to set the number of jobs that should be executed in parallel when running the pipeline. 

This pipeline will download *de novo* germline mutation data for the BXDs, annotate it with relevant metadata, run a genome-wide IHD scan, find any significant markers, and plot the results of the scan.

### Running an inter-haplotype distance scan

Once you have assembled the input files above, a scan can be performed as follows:

```
python scripts/run_ihd_scan.py \
        --singletons /path/to/mutation/csv \
        --config /path/to/config/json \
        --out /name/of/output/csv \
```

There are a small number of optional arguments:

* `-k` sets the kmer size to use for the mutation types (k = 1 will compute distances between aggregate 1-mer mutation spectra, k = 3 will compute distances between aggregate 3-mer mutation spectra). Default value is 1. 

* `-permutations` sets the number of permutations to use when calculating significance thresholds for the IHD scan. Default value is 1,000.

### Plotting the results of an IHD scan

```
python scripts/plot_ihd_results.py \
        --results /path/to/output/csv \
        --markers /path/to/marker/metadata/csv \
        --outpref output_prefix
```

There is one optional argument, `-colname`, that can be used to specify the name of the column in the marker metadata CSV that indicates the physical/genetic map position you wish to plot in the Manhattan plot. The argument defaults to "Mb."

## Running tests

Tests can be run using `pytest` from the root level of the project directory as:

```
pytest .
```

These tests are run automatically via GitHub actions after each commit to this repository, as well. See the badge at the top of the README for the status (pass/fail) of those tests.

## Project layout

    ihd/                                # code for running the IHD method
        utils.py                        # bulk of utility functions
        run_ihd_scan.py                 # wrapper that calls utilities for computing IHD
        plot_ihd_results.py             # script for plotting Manhattan-esque results
        schema.py                       # pandera schema used to validate input/output dataframes

    scripts/                            # various scripts for analyses presented in paper
        combine_bxd_singletons.py       # combine and reformat the downloaded BXD singleton data
        compare_spectra_at_marker.py    # plot mutation spectra in BXDs that inherited either genotype at either of the two mutator loci
        find_nonsynonymous.py           # find any nonsynonymous mutations that are fixed differences between the parental DBA/2J and C57BL/6J strains
        compute_mutation_spectra.py     # generate heatmaps comparing 3-mer mutation spectra between two collections of mutations
        compare_mgp_spectra.py          # plot mutation spectra in Sanger MGP strains that have either genotype at either of the two mutator loci
        test_for_epistasis.R            # test for evidence of epistasis between mutator loci on C>A mutation counts using generalized linear models
        run_pipeline.smk                # snakemake pipeline to generate manuscript figures

    tests/
        fixtures.py                     # fixtures used by `pytest`
        test_utils.py                   # testing suite for methods in `utils.py`

    data/
        genotypes/                      # directory containing formatted `.geno` files for the BXDs that contain sample genotypes at every tested marker
        json/                           # directory containing JSON configuration files for IHD scans using the BXDs
        mutations/                      # directory containing per-sample *de novo* mutation data in the BXDs
        bam_names_to_metadata.xlsx      # Excel file with metadata about the BXD RILs
