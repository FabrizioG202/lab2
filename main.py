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
confusion_matrix = sklearn.metrics.confusion_matrix(
    combined_test["label"].to_list(),
    [
        1 if max(pswm.get_scores(motif[:90])) >= threshold else 0
        for motif in combined_test["sequence"].to_list()
    ],
)


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

    def _calc_params(name: str, values: list[float]) -> dict[str, int | float | Any]:
        return {
            f"max_{name}": np.max(values),
            f"avg_{name}": np.mean(values),
        }

    with contextlib.redirect_stderr(io.StringIO()):
        analysis = ProteinAnalysis(_pad_sequence(sequence))

        out = {
            **{f"comp_{aa}": v for aa, v in get_composition(sequence[:22]).items()},
            **_calc_params(
                "hydrophobicity",
                analysis.protein_scale(Bio.SeqUtils.ProtParamData.kd, 5),
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

    return out


feature_df = pl.DataFrame(
    [
        {
            "accession": row["accession"],
            "label": row["label"],
            **get_features(row["sequence"]),
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
    X, y, test_size=0.2, random_state=42, stratify=y
)


pipeline = sklearn.pipeline.Pipeline(
    [
        ("scaler", sklearn.preprocessing.StandardScaler()),
        ("svm", sklearn.svm.SVC()),
    ]
)

PARAMS = [
    {
        "svm__kernel": ["rbf"],
        "svm__C": [0.1, 1, 10, 100],
        "svm__gamma": ["scale", 0.001, 0.01, 0.1, 1],
    },
    {
        "svm__kernel": ["linear"],
        "svm__C": [0.1, 1, 10, 100],
    },
]

cv = sklearn.model_selection.StratifiedKFold(n_splits=5, shuffle=True, random_state=42)

best_score = -np.inf
best_params = None
best_model = None

# Manual CV validation so that we can use tqdm to track the progress of the search.
for params in tqdm(sklearn.model_selection.ParameterGrid(PARAMS), desc="SVM search"):
    model = sklearn.base.clone(pipeline)
    model.set_params(**params)

    scores = sklearn.model_selection.cross_val_score(
        model, X_train, y_train, cv=cv, scoring="f1_macro", n_jobs=-1
    )

    mean_score = scores.mean()

    if mean_score > best_score:
        best_score = mean_score
        best_params = params
        best_model = sklearn.base.clone(pipeline).set_params(**params)
best_score

_ = best_model.fit(X_train, y_train)  # ty:ignore[unresolved-attribute]

print(sklearn.metrics.confusion_matrix(y_test, best_model.predict(X_test)))  # ty:ignore[unresolved-attribute]

result = sklearn.inspection.permutation_importance(
    best_model,
    X_test,
    y_test,
    scoring="f1_macro",
    n_repeats=10,
    random_state=42,
    n_jobs=-1,
)

pi_importance_df = pl.DataFrame(
    {
        "feature": X_df.columns,
        "importance_mean": result.importances_mean,
        "importance_std": result.importances_std,
    }
).sort("importance_mean", descending=True)


rf = sklearn.ensemble.RandomForestClassifier(
    n_estimators=500, random_state=42, n_jobs=-1, class_weight="balanced"
)

_ = rf.fit(X_train, y_train)

rf_importance_df = pl.DataFrame(
    {
        "feature": X_df.columns,
        "importance": rf.feature_importances_,
    }
).sort("importance", descending=True)

# plot the pi and rf importnace as side by side horizontal bar plots using plotly
fig = go.Figure(
    data=[
        go.Bar(
            x=pi_importance_df["importance_mean"],
            y=pi_importance_df["feature"],
            orientation="h",
            name="Permutation Importance",
            error_x=dict(type="data", array=pi_importance_df["importance_std"]),
        ),
        go.Bar(
            x=rf_importance_df["importance"],
            y=rf_importance_df["feature"],
            orientation="h",
            name="Random Forest Importance",
        ),
    ]
)

fig.update_layout(
    title="Feature Importance",
    xaxis_title="Importance",
    yaxis_title="Feature",
    yaxis={
        "categoryorder": "total ascending",
        "tickfont": {"size": 10},
        "tickmode": "linear",
    },
    barmode="group",
)

fig.show()


# Neural Network Training Pipeline
from sklearn.neural_network import MLPClassifier

nn = MLPClassifier(
    hidden_layer_sizes=(100, 50),
    max_iter=1000,
    random_state=42,
    early_stopping=True,
    validation_fraction=0.1,
)

_ = nn.fit(X_train, y_train)

nn_predictions = nn.predict(X_test)
nn_confusion_matrix = sklearn.metrics.confusion_matrix(y_test, nn_predictions)

print("Neural Network Confusion Matrix:")
print(nn_confusion_matrix)

# Plot the neural network confusion matrix
plot_confusion_matrix(nn_confusion_matrix).write_image(
    ".imgs/neural_network_confusion_matrix.svg"
)
