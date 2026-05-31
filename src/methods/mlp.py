import numpy as np
from sklearn.neural_network import MLPClassifier


def fit_mlp(
    x_train: np.ndarray,
    y_train: np.ndarray,
):

    model = MLPClassifier(
        hidden_layer_sizes=(100, 50),
        max_iter=1000,
        random_state=42,
        early_stopping=True,
        validation_fraction=0.1,
    )

    _ = model.fit(x_train, y_train)

    return model
