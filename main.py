from sklearn.neural_network import MLPClassifier
from sklearn.feature_extraction.text import TfidfVectorizer
from sklearn.model_selection import train_test_split
from tqdm import tqdm
import Bio.SeqUtils.ProtParamData
from Bio.SeqUtils.ProtParam import ProteinAnalysis
import importlib
import itertools
import re
from dataclasses import dataclass
from typing import Iterable, Iterator, Any, Protocol, cast
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
import src.data_collection; importlib.reload(src.data_collection); clear_output()  # noqa: E703, E702, E402 # fmt: skip
import src.logo_generator; importlib.reload(src.logo_generator); clear_output()  # noqa: E703, E702, E402 # fmt: skip
import src.utils.AdditionalProtParamData; importlib.reload(src.utils.AdditionalProtParamData); clear_output()  # noqa: E703, E702, E402 # fmt: skip
import src.methods.von_heijne; importlib.reload(src.methods.von_heijne);  # noqa: E703, E702, E402 # fmt: skip
import src.processing; importlib.reload(src.processing); clear_output()  # noqa: E703, E702, E402 # fmt: skip
import src.graphics; importlib.reload(src.graphics); clear_output()  # noqa: E703, E702, E402 # fmt: skip
import src.feature_extraction; importlib.reload(src.feature_extraction); clear_output()  # noqa: E703, E702, E402 # fmt: skip
import src.methods.svm; importlib.reload(src.methods.svm); clear_output()  # noqa: E703, E702, E402 # fmt: skip
import src.methods.mlp; importlib.reload(src.methods.mlp); clear_output()  # noqa: E703, E702, E402 # fmt: skip
import src.methods.model_selection; importlib.reload(src.methods.model_selection); clear_output()  # noqa: E703, E702, E402 # fmt: skip
import src.methods.von_heijne; importlib.reload(src.methods.von_heijne); clear_output()  # noqa: E703, E702, E402 # fmt: skip
import src.feature_extraction; importlib.reload(src.feature_extraction); clear_output()  # noqa: E703, E702, E402 # fmt: skip
import src.feature_importance; importlib.reload(src.feature_importance); clear_output()  # noqa: E703, E702, E402 # fmt: skip

# collect data
all_positive, all_negative, positive, negative = src.data_collection.collect_data()

# save metrics images.
src.data_collection.save_exploration_plots(all_positive, all_negative)


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
# We assume this does not throw since we selected for residues with SP cleavage site > 14 and sequence length > 40
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
).savefig("report/.imgs/positive_logo.svg")


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

# ██╗   ██╗  ██████╗  ███╗   ██╗        ██╗  ██╗ ███████╗      ██╗ ███╗   ██╗ ███████╗
# ██║   ██║ ██╔═══██╗ ████╗  ██║        ██║  ██║ ██╔════╝      ██║ ████╗  ██║ ██╔════╝
# ██║   ██║ ██║   ██║ ██╔██╗ ██║ █████╗ ███████║ █████╗        ██║ ██╔██╗ ██║ █████╗
# ╚██╗ ██╔╝ ██║   ██║ ██║╚██╗██║ ╚════╝ ██╔══██║ ██╔══╝   ██   ██║ ██║╚██╗██║ ██╔══╝
#  ╚████╔╝  ╚██████╔╝ ██║ ╚████║        ██║  ██║ ███████╗ ╚█████╔╝ ██║ ╚████║ ███████╗
#   ╚═══╝    ╚═════╝  ╚═╝  ╚═══╝        ╚═╝  ╚═╝ ╚══════╝  ╚════╝  ╚═╝  ╚═══╝ ╚══════╝


# This function creates, trains, and evaluates the von heijne model
# and returns the position-specific weight matrix and the confusion matrix for the model.
pswm, von_heijne_confusion_matrix = src.methods.von_heijne.compute_pswm(
    # 80% training split.
    *src.processing.split_df(combined, fraction=0.8, seed=42)
)

# describe the confusion matrix for the von heijne model
print(von_heijne_confusion_matrix.describe())

# Save the confusion matrix for the von heijne model as an image
von_heijne_confusion_matrix.plot().write_image(
    "report/.imgs/von_heijne_confusion_matrix.svg"
)

# ███████╗ ██╗   ██╗ ███╗   ███╗
# ██╔════╝ ██║   ██║ ████╗ ████║
# ███████╗ ██║   ██║ ██╔████╔██║
# ╚════██║ ╚██╗ ██╔╝ ██║╚██╔╝██║
# ███████║  ╚████╔╝  ██║ ╚═╝ ██║
# ╚══════╝   ╚═══╝   ╚═╝     ╚═╝

# This extracts the features, trains the SVM model and evaluates it,
# returning the feature dataframe, the trained model and the confusion matrix for the model.
fit_result = src.methods.svm.fit_svm_models(combined)

# describe the confusion matrix for the full svm model
print(
    f"""Full SVM model confusion matrix:
    {fit_result.full_confusion_matrix.describe()}

    Reduced SVM model confusion matrix:
    {fit_result.reduced_confusion_matrix.describe()}

    Forest Model confusion matrix:
    {fit_result.random_forest_confusion_matrix.describe()}
    """
)


# Save the confusion matrix for the full svm model as an image
fit_result.save_confusion_matrices("report/.imgs")
