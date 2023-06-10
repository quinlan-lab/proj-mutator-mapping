from pandera import Column, Check, DataFrameSchema


IHDResultSchema = DataFrameSchema({
    "marker": Column(str),
    # technically, cosine distance can only take a value from 0 to 2.
    # since we adjust cosine distances and potentially even use chisquare
    # statistics, we won't apply a strict check to the range of possible
    # distance values.
    "Distance": Column(float), 
    "k": Column(int, Check.isin([1, 3])),
})

MutationSchema = DataFrameSchema({
    "sample": Column(str),
    "kmer": Column(str, Check.str_matches(r"[ATCGN]{3}>[ATCGN]{3}")),
    "count": Column(int, Check.greater_than_or_equal_to(0)),
})

MarkerMetadataSchema = DataFrameSchema({
    "marker": Column(str),
    "chromosome": Column(str, required=False, coerce=True),
    "cM": Column(float, required=False),
    "Mb": Column(float, required=False),
})
