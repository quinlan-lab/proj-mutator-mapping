import numpy as np
from scripts.utils import (
    compute_haplotype_distance,
    compute_mean,
    compute_spectra,
)
import pytest

def test_compute_haplotype_distance(
    wt_haplotype_array: np.ndarray,
    mut_haplotype_array: np.ndarray,
):

    
    exp = 0.01872

    assert np.isclose(
        compute_haplotype_distance(
            wt_haplotype_array,
            mut_haplotype_array,
        ),
        exp,
        rtol=1e-4
    )


def test_compute_mean(
    wt_haplotype_array: np.ndarray,
    mut_haplotype_array: np.ndarray,
):
    assert np.all(np.isclose(
        compute_mean(wt_haplotype_array),
        np.array([8.1, 9.7, 11, 7.5, 12.7, 8.4]),
    ))
    assert np.all(np.isclose(
        compute_mean(mut_haplotype_array),
        np.array([8.7, 7.7, 11.3, 11.4, 8, 10.2]),
    ))

@pytest.mark.xfail
def test_compute_spectra(bad_mutation_dataframe):
    assert np.array_equal(
        compute_spectra(bad_mutation_dataframe, k=1),
        np.array([
            [0, 3, 4, 2, 1, 1],
            [0, 3, 4, 2, 1, 1],
            [1, 1, 1, 1, 1, 1],
        ]),
    )

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
