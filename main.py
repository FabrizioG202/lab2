import importlib
import re

import logomaker
import numpy as np
import polars as pl
from cycler import L
from fontTools.merge.layout import first

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


# Setup the working director for the project
# (creates .data, and .imgs)
collector.clean_cache()
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

# K_RESIDUES_BEFORE = 13
# K_RESIDUES_AFTER = 2

K_RESIDUES_BEFORE = 6
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

importlib.reload(src.logo_generator)

# generate the logo
fig = src.logo_generator.generate_logo(
    positive["motif"].to_list(),
    left_bound=K_RESIDUES_BEFORE,
    right_bound=K_RESIDUES_AFTER,
)
fig.savefig(".imgs/positive_logo.svg")

# possible aas
# ALPHABET = "GAVPLIMFWYSTCNQHDEKR"
# ALPHABET = sorted(list("GAVPLIMFWYSTCNQHDEKR"))
ALPHABET = list("ARNDCQEGHILKMFPSTWYV")

# split 80% train 20% test
# positive_train = positive.sample(fraction=0.8, seed=42)

# Create PSWM (Position Score Weight Matrix) as a 2D array with shape (len(ALPHABET), K_RESIDUES_BEFORE + K_RESIDUES_AFTER)
pswm = np.ones((K_RESIDUES_BEFORE + K_RESIDUES_AFTER, len(ALPHABET)), dtype=float)

# motifs = positive_train["motif"].to_list()
motifs = ["STAAQAEP", "AVESSPIF", "LTVALAAE", "LSLSQSTN", "MIGVESVR", "SKPTRAFS"]

# loop over motifs
for motif in motifs:
    for i, aa in enumerate(motif):
        if aa in ALPHABET:
            # TODO: This is really inefficient, we should create a map instead.
            aa_index = ALPHABET.index(aa)
            pswm[i, aa_index] += 1


PWSM = np.log2(
    np.round(
        # Normalize by the number of sequences and 20 pseudocounts
        # (pswm / (len(motifs) + 20))
        np.round(pswm / (len(motifs) + 20), 2)
        / [
            src.utils.AdditionalProtParamData.swissprot_composition_2[aa]
            for aa in ALPHABET
        ],
        1,
    )
)

np.round(PWSM, 1)


def get_scores_for_sequence(
    PWSM: np.ndarray, sequence: str, alphabet: list[str]
) -> list[float]:

    scores = []

    window_size = PWSM.shape[0]

    # We select only the first 90 aa
    seq = sequence[:90]

    for i in range(len(seq) - window_size + 1):
        window = seq[i : i + window_size]
        score = 0

        for j, aa in enumerate(window):
            if aa in alphabet:
                aa_index = alphabet.index(aa)
                score += PWSM[j, aa_index]

        scores.append(score)

    return scores


# TEST_SEQ = positive[0]["sequence"].to_list()[0]

scores = get_scores_for_sequence(PWSM, "MRFLAATFLLLALSTAAQAEPVQF", ALPHABET)
np.max(np.round(scores, 1))
