import pytest
import numpy as np
import pandas as pd

rng = np.random.default_rng(seed = 42)

@pytest.fixture
def wt_haplotype_array() -> np.ndarray:
    """generate a random array of simulated 
    1-mer mutation counts for 10 haplotypes

    Returns:
        np.ndarray
    """
    return rng.integers(low=0, high=20, size=(10, 6)).astype(np.float32)

@pytest.fixture
def mut_haplotype_array() -> np.ndarray:
    """generate a random array of simulated 
    1-mer mutation counts for 10 haplotypes

    Returns:
        np.ndarray
    """
    return rng.integers(low=0, high=20, size=(10, 6)).astype(np.float32)

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
        'Strain': np.repeat(['A', 'B', 'C'], repeats=50),
        'kmer': rng.choice(kmer_mutations, size=150, replace=True),
        'count': [1] * 150,
    })
