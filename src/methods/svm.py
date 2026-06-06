from dataclasses import dataclass

import polars as pl
import sklearn.metrics
from sklearn.base import BaseEstimator

import src.feature_extraction
import src.feature_importance
import src.graphics
import src.methods.model_selection


def _get_svm_models() -> list[sklearn.svm.SVC]:
    """Returns a list of SVM models with different hyperparameters to be used in the grid search."""

    pipeline = sklearn.pipeline.Pipeline(
        [
            ("scaler", sklearn.preprocessing.StandardScaler()),
            ("svm", sklearn.svm.SVC()),
        ]
    )

    PARAMS = [
        {
            "svm__kernel": ["rbf"],
            "svm__C": [0.1, 1, 2, 4, 8],
            "svm__gamma": ["scale", 0.001, 0.01, 0.1, 1, 2],
        },
        {
            "svm__kernel": ["linear"],
            "svm__C": [0.1, 1, 2, 4, 8],
        },
    ]
    models = []

    # Manual CV validation so that we can use tqdm to track the progress of the search.
    for params in sklearn.model_selection.ParameterGrid(PARAMS):
        model = sklearn.base.clone(pipeline)

        model.set_params(**params)
        models.append(sklearn.base.clone(model))

    return models


@dataclass
class SVMFitResult:
    # models
    full_model: sklearn.pipeline.Pipeline
    reduced_model: sklearn.pipeline.Pipeline
    random_forest_model: sklearn.ensemble.RandomForestClassifier

    # confusion matrices
    full_confusion_matrix: src.graphics.ConfusionMatrix
    reduced_confusion_matrix: src.graphics.ConfusionMatrix
    random_forest_confusion_matrix: src.graphics.ConfusionMatrix

    # feature importances
    permutation_importance: pl.DataFrame
    random_forest_importance: pl.DataFrame

    # saves the confustion matrices to the fiven folder
    def save_confusion_matrices(self, folder: str) -> None:
        self.full_confusion_matrix.plot().write_image(
            f"{folder}/full_svm_confusion_matrix.svg"
        )
        self.reduced_confusion_matrix.plot().write_image(
            f"{folder}/reduced_svm_confusion_matrix.svg"
        )
        self.random_forest_confusion_matrix.plot().write_image(
            f"{folder}/random_forest_confusion_matrix.svg"
        )


def fit_svm_models(
    combined: pl.DataFrame, recompute_best_model: bool = False
) -> SVMFitResult:
    # extract the features for the list of sequence
    feature_df = pl.DataFrame(
        [
            {
                "accession": row["accession"],
                "label": row["label"],
                **src.feature_extraction.FeatureExtractor(
                    row["sequence"]
                ).get_all_features(k=22, n=90, window_size=5),
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

    # computed maximizing f1 score we found the best model
    # to run the optimizing process again, uncomment the line below
    best_model = None
    if recompute_best_model:
        best_model, _ = src.methods.model_selection.find_model_by_cv(
            X_train,
            y_train,
            _get_svm_models(),
        )
    else:
        # This is the model we found to be the best after running the optimization process,
        # we hardcode it here to avoid running the optimization process again.
        best_model = sklearn.pipeline.Pipeline(
            [
                ("scaler", sklearn.preprocessing.StandardScaler()),
                (
                    "svm",
                    sklearn.svm.SVC(kernel="rbf", C=2, gamma=0.01, random_state=42),
                ),
            ]
        )

    # refit on the whole training set
    _ = best_model.fit(X_train, y_train)  # ty:ignore[unresolved-attribute]

    # compute the confusion matrix for the best model on the test set and save it as an image
    svm_confusion_matrix = src.graphics.ConfusionMatrix(
        sklearn.metrics.confusion_matrix(
            y_test,
            best_model.predict(X_test),  # ty:ignore[unresolved-attribute]
        )
    )

    # Feature importance analysis using random forest importance and permutation importance on the SVM model
    feature_analyzer = src.feature_importance.FeatureImportance(
        feature_df=feature_df,
        x_train=X_train,
        y_train=y_train,
        best_model=best_model,
    )

    # compute importances using both methods and save the top 5 features for each method as an image
    permutation_importance = feature_analyzer.compute_permutation_importance(
        n_repeats=5
    )
    random_forest_importance, random_forest_model = (
        feature_analyzer.compute_random_forest_importance()
    )

    # Intermezzo, we get the metrics and confusion matrix for the random forest model and save it as an image
    random_forest_confusion_matrix = src.graphics.ConfusionMatrix(
        sklearn.metrics.confusion_matrix(
            y_test,
            random_forest_model.predict(X_test),  # ty:ignore[unresolved-attribute]
        )
    )

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

    # To run the optimizing process again, uncomment the lines below
    reduced_best_model = None
    if recompute_best_model:
        reduced_best_model, _ = src.methods.model_selection.find_model_by_cv(
            X_reduced_train,
            y_reduced_train,
            _get_svm_models(),
        )
    else:
        # This is the model we found to be the best after running the optimization process on the reduced feature set,
        reduced_best_model = sklearn.pipeline.Pipeline(
            [
                ("scaler", sklearn.preprocessing.StandardScaler()),
                ("svm", sklearn.svm.SVC(kernel="rbf", C=2, gamma=0.1, random_state=42)),
            ]
        )

    # fit the model on the reduced training set
    _ = reduced_best_model.fit(X_reduced_train, y_reduced_train)  # ty:ignore[unresolved-attribute]

    # compute and save the confusion matrix for the best model on the reduced feature set on the test set as an image
    reduced_svm_confusion_matrix = src.graphics.ConfusionMatrix(
        sklearn.metrics.confusion_matrix(
            y_reduced_test,
            reduced_best_model.predict(X_reduced_test),  # ty:ignore[unresolved-attribute]
        )
    )

    return SVMFitResult(
        full_model=best_model,
        reduced_model=reduced_best_model,
        random_forest_model=random_forest_model,
        full_confusion_matrix=svm_confusion_matrix,
        reduced_confusion_matrix=reduced_svm_confusion_matrix,
        random_forest_confusion_matrix=random_forest_confusion_matrix,
        permutation_importance=permutation_importance,
        random_forest_importance=random_forest_importance,
    )
