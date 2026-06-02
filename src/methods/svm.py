import contextlib
import importlib
import io
import itertools
import re
from dataclasses import dataclass
from typing import Any, Iterable, Iterator

import Bio.SeqUtils.ProtParamData
import logomaker
import numpy as np
import plotly.express as px
import plotly.graph_objects as go
import plotly.io as pio
import polars as pl
import sklearn.metrics
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from tqdm import tqdm


def get_svm_models(x_train: pl.DataFrame, y_train: pl.Series) -> list[sklearn.svm.SVC]:
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

    # cv = sklearn.model_selection.StratifiedKFold(
    #     n_splits=5, shuffle=True, random_state=42
    # )

    # best_score = -np.inf
    # best_model = None
    models = []

    # Manual CV validation so that we can use tqdm to track the progress of the search.
    for params in tqdm(
        sklearn.model_selection.ParameterGrid(PARAMS), desc="SVM search"
    ):
        model = sklearn.base.clone(pipeline)
        model.set_params(**params)

        models.append(model)
    return models
    # scores = sklearn.model_selection.cross_val_score(
    #     model, x_train, y_train, cv=cv, scoring="f1_macro", n_jobs=-1
    # )

    # mean_score = scores.mean()

    # if mean_score > best_score:
    #     best_score = mean_score
    #     best_model = sklearn.base.clone(pipeline).set_params(**params)

    # _ = best_model.fit(x_train, y_train)  # ty:ignore[unresolved-attribute]

    # return best_model, best_score
    #
