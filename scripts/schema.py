import pandera as pa
from pandera import Column, Check, DataFrameSchema
from typing import Optional


IHDResultSchema = DataFrameSchema({
    "marker": Column(str),
    "distance": Column(float, Check.in_range(min_value=-1, max_value=1)),
    "k": Column(int, Check.isin([1, 3])),
})


MutationSchema = DataFrameSchema({
    "sample": Column(str),
    "kmer": Column(str, Check.str_matches(r"[ATCGN]{3}>[ATCGN]{3}")),
    "count": Column(int, Check.greater_than_or_equal_to(1)),
})


MarkerMetadataSchema = DataFrameSchema({
    "marker": Column(str),
    "cM": Column(float, required=False),
    "Mb": Column(float, required=False),
})
