from typing import Iterable, Union

import numpy as np
import sklearn.base
from tqdm import tqdm


def find_model_by_cv(
    x_train: np.ndarray,
    y_train: np.ndarray,
    candidate_models: Iterable[sklearn.base.BaseEstimator],
    cv: int = 5,
    scoring: str = "f1",
) -> tuple[sklearn.base.BaseEstimator, float]:

    best_model = None
    best_score = -1.0
    models_to_evaluate = list(candidate_models)

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
