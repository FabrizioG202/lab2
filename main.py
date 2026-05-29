import importlib
import itertools
import re
from dataclasses import dataclass
from typing import Iterable, Iterator

import logomaker
import numpy as np
import polars as pl
import sklearn.metrics

import src.data_collection
import src.logo_generator
import src.utils.AdditionalProtParamData

importlib.reload(src.utils.AdditionalProtParamData)
importlib.reload(src.data_collection)

collector = src.data_collection.DataCollector(
    positive_query="""
    (
        (taxonomy_id:2759)
        AND (reviewed:true)
        AND (ft_signal_exp:*)
        AND (fragment:false)
        AND (length:[40 TO *])
        AND (existence:1)
    )
    """,
    negative_query="""
    (
        (reviewed:true)
        AND (fragment:false)
        AND (taxonomy_id:2759)
        AND (length:[40 TO *])
        AND (existence:1)
        AND NOT (ft_signal:*)
        AND (
            (cc_scl_term_exp:SL-0191)
            OR (cc_scl_term_exp:SL-0204)
            OR (cc_scl_term_exp:SL-0039)
            OR (cc_scl_term_exp:SL-0091)
            OR (cc_scl_term_exp:SL-0209)
            OR (cc_scl_term_exp:SL-0173)
        )
    )
    """,
)


collector.clean_cache()

# Setup the working director for the project
# (creates .data, and .imgs)
collector.setup_wd()


# Get the positive examples and cluster them
all_positive = collector.cluster_df(collector.get_positive_examples())

# Get the negative examples and cluster them
all_negative = collector.cluster_df(collector.get_negative_examples())


# positive are only the ones where accession matches cluster_id, so we have one representative per cluster
positive = all_positive.filter(pl.col("accession") == pl.col("cluster_id"))
negative = all_negative.filter(pl.col("accession") == pl.col("cluster_id"))

# Add a column with the sequence neighbouring the Cleavage site.
# ┌───────────┬───────────────────────┬──────────────────────────┐
# │           │                       │                          │
# │           ▼                       ▼                          │
# │LPNTGRLAGCTVFITGASRGIGKAIALKAAKDGANIVIAAKTAQPHPKLLGTIYTAAEEIEA│
# │           ─────────────────────┬───                          │
# │                                ▲                             │
# └────────────────────────────────┼─────────────────────────────┘
#                           cleavage site

K_RESIDUES_BEFORE = 13
K_RESIDUES_AFTER = 2

# Add a column to the positive which has the motifs (the K residues before and after the cleavage site)
# We assume this does not throw since we selected for residues with SP cleavage site > 14 and sequence length > 90
positive = positive.with_columns(
    pl.concat_str(
        [
            pl.col("sequence").str.slice(
                pl.col("cleavage_site") - K_RESIDUES_BEFORE, K_RESIDUES_BEFORE
            ),
            pl.col("sequence").str.slice(pl.col("cleavage_site"), K_RESIDUES_AFTER),
        ]
    ).alias("motif")
)


# generate the logo
fig = src.logo_generator.generate_logo(
    positive["motif"].to_list(),
    left_bound=K_RESIDUES_BEFORE,
    right_bound=K_RESIDUES_AFTER,
)
fig.savefig(".imgs/positive_logo.svg")

# possible aas
ALPHABET = list("GAVPLIMFWYSTCNQHDEKR")

# split 80% train 20% test
positive_train = positive.sample(fraction=0.8, seed=42)


import src.von_heijne; importlib.reload(src.von_heijne);  # noqa: E703 # fmt: skip

pswm = src.von_heijne.PSWM.compute_for(
    positive_train["motif"].to_list(), alphabet=ALPHABET
)

threshold, history = pswm.compute_optimal_threshold(
    folds=5,
    all_sequences=pl.concat(
        [
            positive.with_columns(pl.lit(1).alias("label")).select(
                ["accession", "sequence", "label"]
            ),
            negative.with_columns(pl.lit(0).alias("label")).select(
                ["accession", "sequence", "label"]
            ),
        ]
    ),
)
threshold
