import numpy as np
import plotly.express as px
import plotly.graph_objects as go
import plotly.io as pio


class ConfusionMatrix:
    def __init__(self, confusion_matrix: np.ndarray):
        self.confusion_matrix = confusion_matrix

    @staticmethod
    def from_values(y_true: np.ndarray, y_pred: np.ndarray) -> "ConfusionMatrix":
        from sklearn.metrics import confusion_matrix

        cm = confusion_matrix(y_true, y_pred)
        return ConfusionMatrix(cm)

    @property
    def accuracy(self) -> float:
        tn, fp, fn, tp = self.confusion_matrix.ravel()
        return (tp + tn) / (tp + tn + fp + fn)

    @property
    def precision(self) -> float:
        tn, fp, fn, tp = self.confusion_matrix.ravel()
        return tp / (tp + fp) if (tp + fp) else 0.0

    @property
    def recall(self) -> float:
        tn, fp, fn, tp = self.confusion_matrix.ravel()
        return tp / (tp + fn) if (tp + fn) else 0.0

    @property
    def f1(self) -> float:
        precision = self.precision
        recall = self.recall
        return (
            2 * (precision * recall) / (precision + recall)
            if (precision + recall)
            else 0.0
        )

    @property
    def mcc(self) -> float:
        tn, fp, fn, tp = self.confusion_matrix.ravel()
        denominator = np.sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
        return ((tp * tn) - (fp * fn)) / denominator if denominator else 0.0

    def plot(self) -> "go.Figure":
        return plot_confusion_matrix(self.confusion_matrix)

    def describe(self) -> str:
        return (
            f"Accuracy: {self.accuracy:.4f}\n"
            f"Precision: {self.precision:.4f}\n"
            f"Recall: {self.recall:.4f}\n"
            f"F1 Score: {self.f1:.4f}\n"
            f"MCC: {self.mcc:.4f}"
        )


# plot confusion matrix using plotly
def plot_confusion_matrix(confusion_matrix_data: np.ndarray) -> "go.Figure":
    pio.renderers.default = "png"

    fig = px.imshow(
        confusion_matrix_data,
        labels={"x": "Predicted", "y": "True"},
        x=["Negative", "Positive"],
        y=["Negative", "Positive"],
        text_auto=True,
        color_continuous_scale="sunsetdark",
    )

    fig.update_layout(
        coloraxis_showscale=False,
        height=800,
    )

    return fig
