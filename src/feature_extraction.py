import contextlib
import io
import itertools
from typing import Any

import Bio.SeqUtils.ProtParamData
import numpy as np
import plotly.express as px
import plotly.graph_objects as go
import plotly.io as pio
import sklearn.model_selection
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from scipy.constants import g


class FeatureExtractor:
    sequence: str

    def __init__(self, sequence: str):
        self.sequence = sequence

    @staticmethod
    def _get_aa_composition(sequence: str) -> dict:
        analysis = ProteinAnalysis(sequence)

        return {a: c / len(sequence) for a, c in analysis.count_amino_acids().items()}

    def get_n_term_composition(self, cutoff: int) -> dict:
        """Returns the aa composition up to the n_term_cutoff position in the sequence."""
        return {
            f"n_terminal_comp_{aa}": v
            for aa, v in self._get_aa_composition(self.sequence[:cutoff]).items()
        }

    def get_c_term_composition(self, cutoff: int) -> dict:
        """Returns the aa composition after the n_term_cutoff position in the sequence."""
        return {
            f"c_terminal_comp{aa}": v
            for aa, v in self._get_aa_composition(self.sequence[cutoff:]).items()
        }

    def get_feature_described(
        self, scale: dict[str, float], window_size: int, scale_name: str
    ) -> dict[str, int | float | Any]:

        values = self.get_feature(scale, window_size)
        return {
            f"{scale_name}_mean": np.mean(values),
            f"{scale_name}_std": np.std(values),
            f"{scale_name}_min": np.min(values),
            f"{scale_name}_max": np.max(values),
            f"{scale_name}_max_pos": int(np.argmax(values)),
        }

    def get_feature(
        self,
        scale: dict[str, float],
        window_size: int,
    ) -> list[float]:
        # This is a bit of an hack to prevent noise in the terminal
        # since protein_scale complains about non-standard aminoacids.
        with contextlib.redirect_stderr(io.StringIO()):
            values = ProteinAnalysis(self.sequence).protein_scale(
                scale,
                window=window_size,
            )

        return values


def _generate_feature_extraction_image():
    example_sequence = "MKTAYIAKQRQISFVKSHFSRQDILDLWIYHTQGYFPDWQNYTPGPGIRYPLKFRT"
    extractor = FeatureExtractor(example_sequence)

    values = extractor.get_feature(Bio.SeqUtils.ProtParamData.kd, window_size=5)
    described_features = extractor.get_feature_described(
        Bio.SeqUtils.ProtParamData.kd, window_size=5, scale_name="kd"
    )
    pio.renderers.default = "png"

    # plot values as a line plot
    fig = go.Figure()
    fig.add_trace(
        go.Scatter(
            x=list(range(len(values))),
            y=values,
            mode="lines",
            name="Kyte-Doolittle hydrophobicity",
        )
    )

    # Add horizontal lines for mean, min, max
    fig.add_hline(
        y=described_features["kd_mean"],
        line_dash="dash",
        line_color="green",
        annotation_text="Mean",
    )
    fig.add_hline(
        y=described_features["kd_min"],
        line_dash="dash",
        line_color="blue",
        annotation_text="Min",
    )
    fig.add_hline(
        y=described_features["kd_max"],
        line_dash="dash",
        line_color="red",
        annotation_text="Max",
    )

    # Add vertical line for max position
    fig.add_vline(
        x=described_features["kd_max_pos"],
        line_dash="dash",
        line_color="orange",
        annotation_text="Max Pos",
    )
    fig.add_vline(
        x=22,
        line_color="yellow",
        annotation_text="Average Cut Position",
    )

    fig.update_layout(
        xaxis_title="Windowed Position in sequence",
        yaxis_title="Hydrophobicity",
    )
    fig.show()
    fig.write_image("../report/.imgs/feature_extraction_example.svg")


# _generate_feature_extraction_image()
