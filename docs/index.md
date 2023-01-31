# Mapping mutator alleles with inter-haplotype distance

## Usage

Identify alleles that affect the mutation spectrum in bi-parental recombinant inbred crosses. 

## Project layout

    scripts/
        compute_distance.py  # wrapper that calls utilities for computing inter-haplotype distances (IHD)
        plot_manhattan.py    # code used to plot results of IHD scans
        utils.py             # bulk of the actual methods used for IHD
    tests/
        fixtures.py          # fixtures used by `pytest`
        test_utils.py        # testing suite for methods in `utils.py`
    data/
        genotypes/           # directory containing formatted `.geno` files that contain sample genotypes at every tested marker
        json/                # directory containing JSON configuration files for IHD scans
        singletons/          # directory containing per-sample *de novo* mutation data

