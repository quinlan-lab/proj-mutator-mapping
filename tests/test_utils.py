import numpy as np
from ihd.utils import (
    compute_haplotype_distance,
    perform_ihd_scan,
    perform_permutation_test,
    compute_colmean,
    compute_spectra,
    shuffle_spectra,
    compute_nanrowsum,
    compute_allele_frequency,
    compute_genotype_similarity,
    compute_residuals,
)
import pytest


def test_shuffle_spectra(
    wt_haplotype_array: np.ndarray,
    mut_haplotype_array: np.ndarray,
):
    combined_spectra = np.concatenate(
        (wt_haplotype_array, mut_haplotype_array), )
    shuffled_spectra = shuffle_spectra(combined_spectra)
    assert np.array_equal(combined_spectra, shuffled_spectra) is False


def test_compute_haplotype_distance(
    wt_haplotype_array: np.ndarray,
    mut_haplotype_array: np.ndarray,
):

    exp = 0.08372917

    assert np.isclose(
        compute_haplotype_distance(
            wt_haplotype_array,
            mut_haplotype_array,
        ),
        exp,
    )


def test_compute_mean(
    wt_haplotype_array: np.ndarray,
    mut_haplotype_array: np.ndarray,
):
    assert np.all(np.isclose(
        compute_colmean(wt_haplotype_array),
        np.array([2., 3., 3., 3., 2., 1.]),
    ))
    assert np.all(np.isclose(
        compute_colmean(mut_haplotype_array),
        np.array([2., 2., 5., 3., 4., 3.]),
    ))

@pytest.mark.xfail
def test_compute_spectra(bad_mutation_dataframe):
    _, _, _ = compute_spectra(bad_mutation_dataframe, k=1)


def test_compute_rowsum(genotype_array):
    rowsums = compute_nanrowsum(genotype_array)
    assert np.array_equal(np.array([4., 3., 2.]), rowsums,)

def test_compute_rowsum_nans(genotype_array_nans):
    rowsums = compute_nanrowsum(genotype_array_nans)
    assert np.array_equal(np.array([4., 2., 0.]), rowsums,)

def test_compute_allele_frequency_nans(genotype_array_nans):
    af = compute_allele_frequency(genotype_array_nans)
    assert np.array_equal(np.array([0.5, 0.25, 0.]), af,)


def test_compute_genotype_similarity(genotype_array):
    genotype_sims = compute_genotype_similarity(genotype_array)
    exp = np.array([
        np.corrcoef([1., 0.25, 0.5], [0., 0.5, 0.]),
        np.corrcoef([0., 1., 0.], [0.5, 0., 0.5]),
        np.corrcoef([1., 0., 1.], [1/3, 0.5, 0.]),
    ])


def test_compute_spectra(good_mutation_dataframe):
    samples, mutations, spectra = compute_spectra(good_mutation_dataframe, k=1)
    assert samples == ["A", "B", "C"]
    assert mutations == ["A>C", "A>G", "A>T", "C>A", "C>G", "C>T"]
    assert np.array_equal(
        spectra,
        np.array([
            [11, 9, 7, 14, 9, 0],
            [9, 11, 7, 10, 13, 0],
            [11, 10, 4, 7, 14, 4],
        ]),
    )

@pytest.mark.parametrize("X,y", [(np.arange(5).astype(np.float64), np.arange(5).astype(np.float64)),])
def test_compute_residuals(X, y):
    assert np.allclose(compute_residuals(X, y), np.zeros(5))

def test_perform_ihd_scan(spectra_array, genotype_array, genotype_similarity):
    focal_dists = perform_ihd_scan(
        spectra_array,
        genotype_array,
        genotype_similarity,
    )

    raw_dists = np.array([
                           0.083842,
                           0.273145,
                           0.132533,
                       ])
    assert np.allclose(np.array(focal_dists),
                       raw_dists - np.mean(raw_dists))


@pytest.mark.parametrize("n_permutations,comparison_wide", [
    (100, False),
    (1000, False),
    (100, True),
])
def test_perform_ihd_permutation_test(spectra_array, genotype_array, genotype_similarity, n_permutations, comparison_wide):
    perm_res = perform_permutation_test(
        spectra_array,
        genotype_array,
        genotype_similarity,
        n_permutations=n_permutations,
        comparison_wide=comparison_wide,
    )
    if comparison_wide:
        assert perm_res.shape[0] == n_permutations and perm_res.shape[1] == genotype_array.shape[0]
    else:
        assert perm_res.shape[0] == n_permutations
