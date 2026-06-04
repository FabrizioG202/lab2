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
import src.feature_extraction; importlib.reload(src.feature_extraction); clear_output()  # noqa: E703, E702, E402 # fmt: skip
import src.feature_importance; importlib.reload(src.feature_importance); clear_output()  # noqa: E703, E702, E402 # fmt: skip

# collect data
all_positive, all_negative, positive, negative = src.data_collection.collect_data()

src.data_collection.plot_lengths(
    positive,
    negative,
).write_image("report/.imgs/length_distribution.svg")

src.data_collection.plot_kingdom_distribution(
    positive,
    negative,
).write_image("report/.imgs/kingdom_distribution.svg")

src.data_collection.plot_cleaved_region_lengths(
    positive,
).write_image("report/.imgs/cleaved_region_length_distribution.svg")


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


# split 80% train 20% test
combined_train, combined_test = src.processing.split_df(combined, fraction=0.8, seed=42)

# Compute the PSWM on the training set, using only the positive examples
# with a motif (i.e. where we could extract the motif).
pswm = src.methods.von_heijne.PSWM.compute_for(
    [m for m in combined_train["motif"].to_list() if m is not None],
    alphabet=list("GAVPLIMFWYSTCNQHDEKR"),
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


# Now, evaluate on the test set and compute the confusion matrix for the PSWM method and save it as an image
von_heijne_confusion_matrix = src.graphics.ConfusionMatirx(
    sklearn.metrics.confusion_matrix(
        combined_test["label"].to_list(),
        [
            1 if max(pswm.get_scores(motif[:90])) >= threshold else 0
            for motif in combined_test["sequence"].to_list()
        ],
    )
)

von_heijne_confusion_matrix.plot().write_image(
    "report/.imgs/von_heijne_confusion_matrix.svg"
)

print(von_heijne_confusion_matrix.describe())

# ███████╗ ██╗   ██╗ ███╗   ███╗
# ██╔════╝ ██║   ██║ ████╗ ████║
# ███████╗ ██║   ██║ ██╔████╔██║
# ╚════██║ ╚██╗ ██╔╝ ██║╚██╔╝██║
# ███████║  ╚████╔╝  ██║ ╚═╝ ██║
# ╚══════╝   ╚═══╝   ╚═╝     ╚═╝


# extract the features for the list of sequence
feature_df = pl.DataFrame(
    [
        {
            "accession": row["accession"],
            "label": row["label"],
            **src.feature_extraction.FeatureExtractor(row["sequence"]).get_all_features(
                k=22, n=90, window_size=5
            ),
        }
        for row in combined.iter_rows(named=True)
    ],
)

# Extract the features
X_df = feature_df.drop(["label", "accession"])
X = X_df.to_numpy()
y = feature_df["label"].to_numpy()

# split 80% train 20% test
X_train, X_test, y_train, y_test = sklearn.model_selection.train_test_split(
    X,
    y,
    test_size=0.2,
    random_state=42,
    stratify=y,
)
svm_models = src.methods.svm.get_svm_models()[::-1]


class ClassifierEstimator(Protocol):
    def fit(self, X, y): ...
    def predict(self, X): ...


def find_model_by_cv(
    x_train: np.ndarray,
    y_train: np.ndarray,
    candidate_models: Iterable[sklearn.base.BaseEstimator] | None = None,
    cv: int = 5,
    scoring: str = "f1",
) -> tuple[sklearn.base.BaseEstimator, float]:

    best_model = None
    best_score = -1.0
    models_to_evaluate = (
        list(candidate_models) if candidate_models is not None else svm_models
    )

    for model in tqdm(models_to_evaluate):
        cv_results = sklearn.model_selection.cross_validate(
            model,
            x_train,
            y_train,
            cv=cv,
            scoring=scoring,
            return_train_score=False,
        )

        if (mean_score := np.mean(cv_results["test_score"])) > best_score:
            best_score = mean_score
            best_model = sklearn.base.clone(model)

    return best_model, best_score  # ty:ignore[invalid-return-type]


# computed maximizing f1 score we found the best model
# to run the optimizing process again, uncomment the line below
# best_model = find_model_by_cv()
# This is the model we found to be the best after running the optimization process,
# we hardcode it here to avoid running the optimization process again.
best_model = sklearn.pipeline.Pipeline(
    [
        ("scaler", sklearn.preprocessing.StandardScaler()),
        ("svm", sklearn.svm.SVC(kernel="rbf", C=2, gamma=0.01, random_state=42)),
    ]
)


# refit on the whole training set
_ = best_model.fit(X_train, y_train)

# compute the confusion matrix for the best model on the test set and save it as an image
svm_confusion_matrix = src.graphics.ConfusionMatirx(
    sklearn.metrics.confusion_matrix(
        y_test,
        best_model.predict(
            X_test,
        ),
    )
)
print(svm_confusion_matrix.describe())

svm_confusion_matrix.plot().write_image("report/.imgs/svm_confusion_matrix.svg")


#  █████╗  ███╗   ██╗  █████╗  ██╗  ██╗   ██╗ ███████╗ ██╗ ███████╗
# ██╔══██╗ ████╗  ██║ ██╔══██╗ ██║  ╚██╗ ██╔╝ ██╔════╝ ██║ ██╔════╝
# ███████║ ██╔██╗ ██║ ███████║ ██║   ╚████╔╝  ███████╗ ██║ ███████╗
# ██╔══██║ ██║╚██╗██║ ██╔══██║ ██║    ╚██╔╝   ╚════██║ ██║ ╚════██║
# ██║  ██║ ██║ ╚████║ ██║  ██║ ███████╗██║    ███████║ ██║ ███████║
# ╚═╝  ╚═╝ ╚═╝  ╚═══╝ ╚═╝  ╚═╝ ╚══════╝╚═╝    ╚══════╝ ╚═╝ ╚══════╝

# Feature importance analysis using random forest importance and permutation importance on the SVM model
feature_analyzer = src.feature_importance.FeatureImportance(
    feature_df=feature_df,
    x_train=X_train,
    y_train=y_train,
    best_model=best_model,
)

# compute importances using both methods and save the top 5 features for each method as an image
permutation_importance = feature_analyzer.compute_permutation_importance(
    n_repeats=5,
)
random_forest_importance, random_forest_model = (
    feature_analyzer.compute_random_forest_importance()
)

print("Top 5 features by permutation importance:")
print(permutation_importance.head(5))
print("Top 5 features by random forest importance:")
print(random_forest_importance.head(5))

# Intermezzo, we get the metrics and confusion matrix for the random forest model and save it as an image
random_forest_confusion_matrix = src.graphics.ConfusionMatirx(
    sklearn.metrics.confusion_matrix(
        y_test,
        random_forest_model.predict(X_test),
    )
)

print(random_forest_confusion_matrix.describe())


# Re-select the best SVM model on the top 5 features from permutation importance
selected_features = set(
    permutation_importance["feature"].head(5).to_list()
    + random_forest_importance["feature"].head(5).to_list()
)
X_reduced_df = X_df.select(selected_features)
X_reduced = X_reduced_df.to_numpy()

X_reduced_train, X_reduced_test, y_reduced_train, y_reduced_test = (
    sklearn.model_selection.train_test_split(
        X_reduced,
        y,
        test_size=0.2,
        random_state=42,
        stratify=y,
    )
)

reduced_best_model, _ = find_model_by_cv(
    X_reduced_train,
    y_reduced_train,
)
_ = reduced_best_model.fit(X_reduced_train, y_reduced_train)  # ty:ignore[unresolved-attribute]

reduced_svm_confusion_matrix = src.graphics.ConfusionMatirx(
    sklearn.metrics.confusion_matrix(
        y_reduced_test,
        reduced_best_model.predict(X_reduced_test),  # ty:ignore[unresolved-attribute]
    )
)
print(reduced_svm_confusion_matrix.describe())
reduced_svm_confusion_matrix.plot().write_image(
    "report/.imgs/reduced_svm_confusion_matrix.svg"
)
