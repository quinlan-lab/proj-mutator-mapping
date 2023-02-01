import pytest
import numpy as np
import pandas as pd

rng = np.random.default_rng(seed = 42)

@pytest.fixture
def wt_haplotype_array() -> np.ndarray:
    """1-mer mutation counts 
    for a toy population of 3 haplotypes

    Returns:
        np.ndarray
    """
    return np.array([
        [1, 2, 2, 4, 2, 2],
        [3, 6, 4, 5, 2, 1],
        [2, 1, 3, 0, 2, 0],
    ]).astype(np.float64)
    

@pytest.fixture
def mut_haplotype_array() -> np.ndarray:
    """1-mer mutation counts 
    for a toy population of 3 haplotypes

    Returns:
        np.ndarray
    """
    return np.array([
        [0, 1, 4, 2, 0, 6],
        [1, 1, 9, 4, 4, 1],
        [5, 4, 2, 3, 8, 2],
    ]).astype(np.float64)

@pytest.fixture
def bad_mutation_dataframe(seed: int = 42) -> pd.DataFrame:

    rng = np.random.default_rng(seed = seed)

    mutations = ["C>A", "C>T", "C>G", "A>T", "A>C", "A>G"]

    return pd.DataFrame({
        'Strain': np.repeat(['A', 'B', 'C'], repeats=8),
        'kmer': rng.choice(mutations, size=24, replace=True),
        'count': [1] * 24,
    })

@pytest.fixture
def good_mutation_dataframe(seed: int = 42) -> pd.DataFrame:

    rng = np.random.default_rng(seed = seed)

    nucs = ["A", "T", "C", "G"]
    mutations = ["C>A", "C>T", "C>G", "A>T", "A>C", "A>G"]
    kmer_mutations = []
    for _ in range(24):
        m = rng.choice(mutations, size=1)[0]
        orig, new = m.split(">")
        fp = rng.choice(nucs, size=1)[0]
        tp = rng.choice(nucs, size=1)[0]
        kmer_mutations.append(f"{fp}{orig}{tp}>{fp}{new}{tp}")


    return pd.DataFrame({
        'sample': np.repeat(['A', 'B', 'C'], repeats=50),
        'kmer': rng.choice(kmer_mutations, size=150, replace=True),
        'count': [1] * 150,
    })
