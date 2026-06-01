from tqdm import tqdm
import Bio.SeqUtils.ProtParamData
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
import contextlib
import io
from IPython.display import clear_output
import src.data_collection
import src.logo_generator
import src.utils.AdditionalProtParamData
import src.methods.von_heijne; importlib.reload(src.methods.von_heijne);  # noqa: E703, E702, E402 # fmt: skip
import src.processing; importlib.reload(src.processing); clear_output()  # noqa: E703, E702, E402 # fmt: skip
import src.graphics; importlib.reload(src.graphics); clear_output()  # noqa: E703, E702, E402 # fmt: skip
import src.feature_extraction; importlib.reload(src.feature_extraction); clear_output()  # noqa: E703, E702, E402 # fmt: skip
import src.methods.svm; importlib.reload(src.methods.svm); clear_output()  # noqa: E703, E702, E402 # fmt: skip
import src.methods.mlp; importlib.reload(src.methods.mlp); clear_output()  # noqa: E703, E702, E402 # fmt: skip
import src.feature_extraction; importlib.reload(src.feature_extraction); clear_output()  # noqa: E703, E702, E402 # fmt: skip
import src.feature_importance; importlib.reload(src.feature_importance); clear_output()  # noqa: E703, E702, E402 # fmt: skip

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

# These contain information about the clusters, but there still are multiple sequences per cluster.
all_positive = collector.cluster_df(collector.get_positive_examples())
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
combined_train, combined_test = src.processing.split_df(combined, fraction=0.8, seed=42)

# ██╗   ██╗  ██████╗  ███╗   ██╗        ██╗  ██╗ ███████╗      ██╗ ███╗   ██╗ ███████╗
# ██║   ██║ ██╔═══██╗ ████╗  ██║        ██║  ██║ ██╔════╝      ██║ ████╗  ██║ ██╔════╝
# ██║   ██║ ██║   ██║ ██╔██╗ ██║ █████╗ ███████║ █████╗        ██║ ██╔██╗ ██║ █████╗
# ╚██╗ ██╔╝ ██║   ██║ ██║╚██╗██║ ╚════╝ ██╔══██║ ██╔══╝   ██   ██║ ██║╚██╗██║ ██╔══╝
#  ╚████╔╝  ╚██████╔╝ ██║ ╚████║        ██║  ██║ ███████╗ ╚█████╔╝ ██║ ╚████║ ███████╗
#   ╚═══╝    ╚═════╝  ╚═╝  ╚═══╝        ╚═╝  ╚═╝ ╚══════╝  ╚════╝  ╚═╝  ╚═══╝ ╚══════╝

pswm = src.methods.von_heijne.PSWM.compute_for(
    [m for m in combined_train["motif"].to_list() if m is not None], alphabet=ALPHABET
)

# Use CV to compute the optimal threshold for the PSWM scores on the training set.
threshold, history = pswm.compute_optimal_threshold(
    folds=5,
    all_sequences=combined_train.select(
        [
            "accession",
            "sequence",
            "label",
        ]
    ),
)

# Now, evaluate on the test set
src.graphics.plot_confusion_matrix(
    sklearn.metrics.confusion_matrix(
        combined_test["label"].to_list(),
        [
            1 if max(pswm.get_scores(motif[:90])) >= threshold else 0
            for motif in combined_test["sequence"].to_list()
        ],
    )
).write_image("report/.imgs/von_heijne_confusion_matrix.svg")


# ███████╗ ███████╗  █████╗  ████████╗ ██╗   ██╗ ██████╗  ███████╗ ███████╗
# ██╔════╝ ██╔════╝ ██╔══██╗ ╚══██╔══╝ ██║   ██║ ██╔══██╗ ██╔════╝ ██╔════╝
# █████╗   █████╗   ███████║    ██║    ██║   ██║ ██████╔╝ █████╗   ███████╗
# ██╔══╝   ██╔══╝   ██╔══██║    ██║    ██║   ██║ ██╔══██╗ ██╔══╝   ╚════██║
# ██║      ███████╗ ██║  ██║    ██║    ╚██████╔╝ ██║  ██║ ███████╗ ███████║
# ╚═╝      ╚══════╝ ╚═╝  ╚═╝    ╚═╝     ╚═════╝  ╚═╝  ╚═╝ ╚══════╝ ╚══════╝

feature_df = pl.DataFrame(
    [
        {
            "accession": row["accession"],
            "label": row["label"],
            **src.feature_extraction.get_all_features(
                row["sequence"],
                alphabet=ALPHABET,
            ),
        }
        for row in tqdm(
            combined.iter_rows(named=True),
            total=combined.height,
            desc="Extracting features",
        )
    ],
)

# SVM Training Pipeline
X_df = feature_df.drop(["label", "accession"])
X = X_df.to_numpy()
y = feature_df["label"].to_numpy()

X_train, X_test, y_train, y_test = sklearn.model_selection.train_test_split(
    X,
    y,
    test_size=0.2,
    random_state=42,
    stratify=y,
)

# ███████╗ ██╗   ██╗ ███╗   ███╗
# ██╔════╝ ██║   ██║ ████╗ ████║
# ███████╗ ██║   ██║ ██╔████╔██║
# ╚════██║ ╚██╗ ██╔╝ ██║╚██╔╝██║
# ███████║  ╚████╔╝  ██║ ╚═╝ ██║
# ╚══════╝   ╚═══╝   ╚═╝     ╚═╝

best_svm_model, _ = src.methods.svm.fit_svm(X_train, y_train)

# Save the confusion matrix computed on the test set as an image
src.graphics.plot_confusion_matrix(
    sklearn.metrics.confusion_matrix(y_test, best_svm_model.predict(X_test))
).write_image("report/.imgs/svm_confusion_matrix.svg")


#  █████╗  ███╗   ██╗  █████╗  ██╗  ██╗   ██╗ ███████╗ ██╗ ███████╗
# ██╔══██╗ ████╗  ██║ ██╔══██╗ ██║  ╚██╗ ██╔╝ ██╔════╝ ██║ ██╔════╝
# ███████║ ██╔██╗ ██║ ███████║ ██║   ╚████╔╝  ███████╗ ██║ ███████╗
# ██╔══██║ ██║╚██╗██║ ██╔══██║ ██║    ╚██╔╝   ╚════██║ ██║ ╚════██║
# ██║  ██║ ██║ ╚████║ ██║  ██║ ███████╗██║    ███████║ ██║ ███████║
# ╚═╝  ╚═╝ ╚═╝  ╚═══╝ ╚═╝  ╚═╝ ╚══════╝╚═╝    ╚══════╝ ╚═╝ ╚══════╝

# Feature importance analysis using random forest importance and permutation importance on the SVM model

# Compute the feature importance for the SVM model and save it as an image
src.feature_importance.draw_feature_importance(
    svm_model=best_svm_model, x_train=X_train, y_train=y_train, x_dataframe=X_df
).write_image("report/.imgs/svm_feature_importance.svg")

# ███╗   ███╗ ██╗      ██████╗
# ████╗ ████║ ██║      ██╔══██╗
# ██╔████╔██║ ██║      ██████╔╝
# ██║╚██╔╝██║ ██║      ██╔═══╝
# ██║ ╚═╝ ██║ ███████╗ ██║
# ╚═╝     ╚═╝ ╚══════╝ ╚═╝

mlp_model = src.methods.mlp.fit_mlp(X_train, y_train)
