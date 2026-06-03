from sklearn.feature_extraction.text import TfidfVectorizer
from sklearn.model_selection import train_test_split
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
import src.data_collection; importlib.reload(src.data_collection); clear_output()  # noqa: E703, E702, E402 # fmt: skip
import src.logo_generator; importlib.reload(src.logo_generator); clear_output()  # noqa: E703, E702, E402 # fmt: skip
import src.utils.AdditionalProtParamData; importlib.reload(src.utils.AdditionalProtParamData); clear_output()  # noqa: E703, E702, E402 # fmt: skip
import src.methods.von_heijne; importlib.reload(src.methods.von_heijne);  # noqa: E703, E702, E402 # fmt: skip
import src.processing; importlib.reload(src.processing); clear_output()  # noqa: E703, E702, E402 # fmt: skip
import src.graphics; importlib.reload(src.graphics); clear_output()  # noqa: E703, E702, E402 # fmt: skip
import src.feature_extraction; importlib.reload(src.feature_extraction); clear_output()  # noqa: E703, E702, E402 # fmt: skip
import src.methods.svm; importlib.reload(src.methods.svm); clear_output()  # noqa: E703, E702, E402 # fmt: skip
import src.methods.mlp; importlib.reload(src.methods.mlp); clear_output()  # noqa: E703, E702, E402 # fmt: skip
import src.feature_extraction; importlib.reload(src.feature_extraction); clear_output()  # noqa: E703, E702, E402 # fmt: skip
import src.feature_importance; importlib.reload(src.feature_importance); clear_output()  # noqa: E703, E702, E402 # fmt: skip

# collect data
all_positive, all_negative, positive, negative = src.data_collection.collect_data()

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
).savefig("report/.imgs/positive_logo.svg")

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

# ██╗   ██╗  ██████╗  ███╗   ██╗        ██╗  ██╗ ███████╗      ██╗ ███╗   ██╗ ███████╗
# ██║   ██║ ██╔═══██╗ ████╗  ██║        ██║  ██║ ██╔════╝      ██║ ████╗  ██║ ██╔════╝
# ██║   ██║ ██║   ██║ ██╔██╗ ██║ █████╗ ███████║ █████╗        ██║ ██╔██╗ ██║ █████╗
# ╚██╗ ██╔╝ ██║   ██║ ██║╚██╗██║ ╚════╝ ██╔══██║ ██╔══╝   ██   ██║ ██║╚██╗██║ ██╔══╝
#  ╚████╔╝  ╚██████╔╝ ██║ ╚████║        ██║  ██║ ███████╗ ╚█████╔╝ ██║ ╚████║ ███████╗
#   ╚═══╝    ╚═════╝  ╚═╝  ╚═══╝        ╚═╝  ╚═╝ ╚══════╝  ╚════╝  ╚═╝  ╚═══╝ ╚══════╝

# split 80% train 20% test
combined_train, combined_test = src.processing.split_df(combined, fraction=0.8, seed=42)

# Compute the PSWM on the training set, using only the positive examples with a motif (i.e. where we could extract the motif).
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
average_cut_site_position = int(
    np.mean(combined_train["cleavage_site"].drop_nulls().to_numpy())
)

best_hp = None
best_f1 = -1.0

for cutoff, model in tqdm(
    src.feature_extraction.get_all_hp_combinations(
        average_cut_position=average_cut_site_position
    )[:5],
    position=0,
):
    feature_df = pl.DataFrame(
        [
            {
                "accession": row["accession"],
                "label": row["label"],
                #
                # composition up to cutoff
                **src.feature_extraction.FeatureExtractor(
                    row["sequence"][:90]
                ).get_c_term_composition(cutoff=cutoff),
                #
                # composition from cutoff
                **src.feature_extraction.FeatureExtractor(
                    row["sequence"][:90]
                ).get_n_term_composition(cutoff=cutoff),
                #
                #
                # global hydrophobicity
                **src.feature_extraction.FeatureExtractor(
                    row["sequence"][:90]
                ).get_feature_described(
                    scale=Bio.SeqUtils.ProtParamData.kd,
                    scale_name="kd",
                    window_size=5,
                ),
                #
                # global alpha-helix tendency
                **src.feature_extraction.FeatureExtractor(
                    row["sequence"][:90]
                ).get_feature_described(
                    scale=src.utils.AdditionalProtParamData.alpha_helix_tendency,
                    scale_name="alpha_helix_tendency",
                    window_size=5,
                ),
                #
                # global transmembrane tendency
                **src.feature_extraction.FeatureExtractor(
                    row["sequence"][:90]
                ).get_feature_described(
                    scale=src.utils.AdditionalProtParamData.transmemberane_tendency,
                    scale_name="transmembrane_tendency",
                    window_size=5,
                ),
                #
                # global bulkiness
                **src.feature_extraction.FeatureExtractor(
                    row["sequence"][:90]
                ).get_feature_described(
                    scale=src.utils.AdditionalProtParamData.bulkiness,
                    scale_name="bulkiness",
                    window_size=5,
                ),
                #
                # global polarity
                **src.feature_extraction.FeatureExtractor(
                    row["sequence"][:90]
                ).get_feature_described(
                    scale=src.utils.AdditionalProtParamData.polarity,
                    scale_name="polarity",
                    window_size=5,
                ),
            }
            for row in combined.iter_rows(named=True)
        ],
    )

    # Extract the features
    X_df = feature_df.drop(["label", "accession"])
    X = X_df.to_numpy()
    y = feature_df["label"].to_numpy()

    best_model: sklearn.base.BaseEstimator = None  # ty:ignore[invalid-assignment]
    best_f1 = -1.0

    # split 80% train 20% test
    X_train, X_test, y_train, y_test = sklearn.model_selection.train_test_split(
        X, y, test_size=0.2, random_state=42, stratify=y
    )

    # CV on the training set to find the best hyperparameters for the model
    cv_results = sklearn.model_selection.cross_validate(
        model,
        X_train,
        y_train,
        cv=5,
        scoring="f1",
        return_train_score=False,
    )

    mean_f1 = np.mean(cv_results["test_score"])

    if mean_f1 > best_f1:
        best_f1 = mean_f1
        best_hp = (cutoff, model)

sequences = combined["sequence"].to_list()
labels = combined["label"].to_list()

X_train, X_test, y_train, y_test = train_test_split(
    sequences, labels, test_size=0.2, stratify=labels, random_state=42
)

vectorizer = TfidfVectorizer(analyzer="char", ngram_range=(3, 6), min_df=2)

X_train_vec = vectorizer.fit_transform(X_train)
X_test_vec = vectorizer.transform(X_test)

model = sklearn.svm.SVC(kernel="rbf", C=1, gamma="scale")
model.fit(X_train_vec, y_train)


# ███████╗ ██╗   ██╗ ███╗   ███╗
# ██╔════╝ ██║   ██║ ████╗ ████║
# ███████╗ ██║   ██║ ██╔████╔██║
# ╚════██║ ╚██╗ ██╔╝ ██║╚██╔╝██║
# ███████║  ╚████╔╝  ██║ ╚═╝ ██║
# ╚══════╝   ╚═══╝   ╚═╝     ╚═╝

# best_svm_model, _ = src.methods.svm.fit_svm(X_train, y_train)

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
