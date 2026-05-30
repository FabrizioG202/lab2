import Bio.SeqUtils.ProtParamData
from curses import window
from Bio.SeqUtils.ProtParam import ProteinAnalysis
import importlib
import itertools
import re
from dataclasses import dataclass
from typing import Iterable, Iterator, Any
import sklearn.metrics
import logomaker
import numpy as np
import polars as pl
import sklearn.metrics
import plotly.express as px
import plotly.graph_objects as go
import plotly.io as pio

import src.data_collection
import src.logo_generator
import src.utils.AdditionalProtParamData
import src.von_heijne; importlib.reload(src.von_heijne);  # noqa: E703, E702 # fmt: skip

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


def split_df(
    df: pl.DataFrame, fraction: float, seed: int
) -> tuple[pl.DataFrame, pl.DataFrame]:
    # create a shuffled version of the dataframe
    shuffled = df.sample(fraction=1.0, shuffle=True, seed=seed)

    # split it into 2
    count = int(len(df) * fraction)
    return shuffled[:count], shuffled[count:]


# combined dataframes with:
# label = 1 for positive
# label = 0 for negative
combined = pl.concat(
    [
        positive.with_columns(pl.lit(1).alias("label")),
        negative.with_columns(pl.lit(0).alias("label")),
    ],
    how="diagonal",
)

# split 80% train 20% test
combined_train, combined_test = split_df(combined, fraction=0.8, seed=42)

# Create the PSWM using the train split.
pswm = src.von_heijne.PSWM.compute_for(
    [m for m in combined_train["motif"].to_list() if m is not None], alphabet=ALPHABET
)

threshold, history = pswm.compute_optimal_threshold(
    folds=5, all_sequences=combined_train.select(["accession", "sequence", "label"])
)
# Now, evaluate on the test set
y_true = combined_test["label"].to_list()
y_pred = [
    1 if max(pswm.get_scores(motif[:90])) >= threshold else 0
    for motif in combined_test["sequence"].to_list()
]


confusion_matrix = sklearn.metrics.confusion_matrix(y_true, y_pred)


# plot confusion matrix using plotly
def plot_confusion_matrix(confusion_matrix: np.ndarray) -> "go.Figure":

    # to make it work with zed's REPL.
    # see: https://github.com/zed-industries/zed/issues/19890
    pio.renderers.default = "png"

    fig = px.imshow(
        confusion_matrix,
        labels={"x": "Predicted", "y": "True"},
        x=["Negative", "Positive"],
        y=["Negative", "Positive"],
        text_auto=True,
        color_continuous_scale="Blues",
    )
    fig.update_layout(title="Confusion Matrix")
    return fig


plot_confusion_matrix(confusion_matrix).write_image(
    ".imgs/von_heijne_confusion_matrix.svg"
)


# A function to get the compositions in aa for the given sequence
# TODO: Make the alphabet a parameter...
def get_composition(sequence: str) -> dict[str, float]:
    composition = {aa: 0 for aa in ALPHABET}
    for aa in sequence:
        if aa in composition:
            composition[aa] += 1
    total = len(sequence)
    return {aa: count / total for aa, count in composition.items()}


def get_features(sequence: str) -> dict[str, float]:
    window_size = 5

    # pad the sequence with X so that with the given
    # window size, we can always get a full window of residues around the cleavage site.
    def _pad_sequence(seq: str) -> str:
        padding = "X" * (window_size // 2)
        return padding + seq + padding

    analysis = ProteinAnalysis(_pad_sequence(sequence))

    def _calc_params(name: str, values: list[float]) -> dict[str, int | float | Any]:
        return {
            f"max_{name}": np.max(values),
            f"avg_{name}": np.mean(values),
        }

    return {
        **{f"comp_{aa}": v for aa, v in get_composition(sequence).items()},
        **_calc_params(
            "hydrophobicity", analysis.protein_scale(Bio.SeqUtils.ProtParamData.kd, 5)
        ),
        **_calc_params(
            "alpha_helix_tendency",
            analysis.protein_scale(
                src.utils.AdditionalProtParamData.alpha_helix_tendency, 5
            ),
        ),
        **_calc_params(
            "transmembrane_tendency",
            analysis.protein_scale(
                src.utils.AdditionalProtParamData.transmemberane_tendency, 5
            ),
        ),
        **_calc_params(
            "bulkiness",
            analysis.protein_scale(src.utils.AdditionalProtParamData.bulkiness, 5),
        ),
        **_calc_params(
            "charge",
            analysis.protein_scale(src.utils.AdditionalProtParamData.polarity, 5),
        ),
    }


TEST_SEQUENCE = positive["sequence"][0]
get_features(TEST_SEQUENCE)
