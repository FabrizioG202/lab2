import numpy as np
import plotly.express as px
import plotly.io as pio


# plot confusion matrix using plotly
def plot_confusion_matrix(confusion_matrix: np.ndarray) -> "go.Figure":
    pio.renderers.default = "png"

    tn, fp, fn, tp = confusion_matrix.ravel()

    accuracy = (tp + tn) / (tp + tn + fp + fn)
    precision = tp / (tp + fp) if (tp + fp) else 0.0
    recall = tp / (tp + fn) if (tp + fn) else 0.0  # sensitivity / true positive rate
    f1 = (
        2 * (precision * recall) / (precision + recall) if (precision + recall) else 0.0
    )

    denominator = np.sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
    mcc = ((tp * tn) - (fp * fn)) / denominator if denominator else 0.0

    fig = px.imshow(
        confusion_matrix,
        labels={"x": "Predicted", "y": "True"},
        x=["Negative", "Positive"],
        y=["Negative", "Positive"],
        text_auto=True,
        color_continuous_scale="sunsetdark",
    )

    metrics_text = (
        f"Accuracy: {accuracy:.3f} | Precision: {precision:.3f} | Recall: {recall:.3f}<br>"
        f"F1: {f1:.3f} | MCC: {mcc:.3f}"
    )

    fig.update_layout(
        coloraxis_showscale=False,
        height=800,
    )

    fig.add_annotation(
        text=metrics_text,
        xref="paper",
        yref="paper",
        x=0.5,
        y=-0.25,
        showarrow=False,
        font=dict(size=14),
        align="center",
    )

    return fig
