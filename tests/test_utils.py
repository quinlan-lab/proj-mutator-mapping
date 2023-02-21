import numpy as np
from ihd.utils import (
    compute_haplotype_distance,
    compute_colmean,
    compute_spectra,
    shuffle_spectra,
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
