import numpy as np
import plotly.express as px
import plotly.io as pio


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
