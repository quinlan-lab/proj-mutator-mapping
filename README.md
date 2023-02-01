# Mapping mutator alleles with inter-haplotype distance

Click the badge to see detailed documentation:

[![Documentation Status](https://readthedocs.org/projects/ansicolortags/badge/?version=latest)](https://quinlan-lab.github.io/proj-mutator-mapping/reference/)

## Summary

Identify alleles that affect the mutation spectrum in bi-parental recombinant inbred crosses. 

## Usage

Before running an inter-haplotype distance (IHD) scan, you'll need to prepare a
small number of input files.

### Description of input files 

1. ***De novo* germline mutation data**

    Mutation data should be in a CSV file with **three required columns** as follows:

    ```
    Strain,kmer,count
    sample_A,CCT>CAT,1
    sample_B,GCA>GAA,1
    ```

    The CSV file can contain either a) one row for every individual mutation observed in each sample, in which case the `count` column should always be 
    set to 1 or b) the aggregate count of every mutation type observed in the sample, in which case the `count` column will reflect the total number of
    each mutation type observed in the sample.

    The dataframe can have any number of additional columns, but only the three 
    defined above will be used.

2. **Marker genotypes**

    Genotype data should be formatted in a similar fashion as in [R/qtl2](https://kbroman.org/qtl2/). Genotypes should be in a CSV file with N rows, where N is the number of genotyped markers. There should be a single column denoting the marker name, and as many columns as there are samples. See below:

    ```
    marker,sample_A,sample_B,sample_C
    rs0001,B,B,D
    rs0002,H,D,B
    rs0003,D,D,D
    ```

3. **Marker information (optional)**

    If you wish to generate Manhattan-esque plots that summarize the results
    of an IHD scan, you'll need to provide a final CSV that links marker IDs with
    either physical or genetic map positions (or both). This file should contain a column called `marker`, and a column specifying one or both of `cM` or `Mb`.

    ```
    marker,cM,Mb
    rs0001,1,4.5230
    rs0002,2.6,5.1994
    rs0003,2.8,5.4872
    ```


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

    The `genotypes` dictionary should map the observed genotypes in file #2 to integer values that will be used during the IHD scan.

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

## Project layout

    scripts/
        run_ihd_scan.py  # wrapper that calls utilities for computing inter-haplotype distances (IHD)
        plot_ihd_scan.py    # code used to plot results of IHD scans
        utils.py             # bulk of the actual methods used for IHD
        schema.py            # pandera schema used to validate dataframes
    tests/
        fixtures.py          # fixtures used by `pytest`
        test_utils.py        # testing suite for methods in `utils.py`
    data/
        genotypes/           # directory containing formatted `.geno` files that contain sample genotypes at every tested marker
        json/                # directory containing JSON configuration files for IHD scans
        mutations/          # directory containing per-sample *de novo* mutation data
