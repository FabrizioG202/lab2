import sklearn.metrics


def get_svm_models() -> list[sklearn.svm.SVC]:
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
            "svm__C": [0.1, 1, 10, 100],
            "svm__gamma": ["scale", 0.001, 0.01, 0.1, 1],
        },
        {
            "svm__kernel": ["linear"],
            "svm__C": [0.1, 1, 10, 100],
        },
    ]
    models = []

    # Manual CV validation so that we can use tqdm to track the progress of the search.
    for params in sklearn.model_selection.ParameterGrid(PARAMS):
        model = sklearn.base.clone(pipeline)

        model.set_params(**params)
        models.append(sklearn.base.clone(model))

    return models
