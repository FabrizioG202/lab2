import contextlib
import io
import itertools
from this import s
from typing import Any

import Bio.SeqUtils.ProtParamData
import numpy as np
import plotly.express as px
import plotly.graph_objects as go
import sklearn.model_selection
from Bio.SeqUtils.ProtParam import ProteinAnalysis


class FeatureExtractor:
    sequence: str

    @staticmethod
    def putative_cutoffs(
        average_cut_position: int,
        step: int = 4,
        first_position: int = 6,
        last_position: int = 40,
    ):
        # similar to a linspace, but ensures that the average_cut_position
        # is included and that the step is consistent
        return list(
            set(
                list(range(first_position, average_cut_position, step))
                + [average_cut_position]
                + list(range(average_cut_position + step, last_position, step))
            )
        )

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

    @staticmethod
    def _get_feature_described(
        sequence: str, scale: dict[str, float], window_size: int, scale_name: str
    ) -> dict[str, int | float]:

        # This is a bit of an hack to prevent noise in the terminal
        # since protein_scale complains about non-standard aminoacids.
        with contextlib.redirect_stderr(io.StringIO()):
            values = ProteinAnalysis(sequence).protein_scale(
                scale,
                window=window_size,
            )
            return {
                f"{scale_name}_mean": np.mean(values),
                f"{scale_name}_std": np.std(values),
                f"{scale_name}_min": np.min(values),
                f"{scale_name}_max": np.max(values),
                f"{scale_name}_max_pos": int(np.argmax(values)),
            }

    def get_feature_described(
        self, scale: dict[str, float], window_size: int, scale_name: str
    ) -> dict[str, int | float]:
        return self._get_feature_described(
            self.sequence,
            scale,
            window_size,
            scale_name,
        )

    def get_n_term_feature_described(
        self, scale: dict[str, float], window_size: int, cutoff: int, scale_name: str
    ) -> dict[str, int | float]:

        return self._get_feature_described(
            self.sequence[:cutoff],
            scale,
            window_size,
            scale_name,
        )

    def get_c_term_feature_described(
        self, scale: dict[str, float], window_size: int, cutoff: int, scale_name: str
    ) -> dict[str, int | float]:

        return self._get_feature_described(
            self.sequence[cutoff:], scale, window_size, scale_name
        )


def get_putative_svm_models() -> list[sklearn.svm.SVC]:
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


def get_all_hp_combinations(
    average_cut_position: int = 22,
) -> list[tuple[int, sklearn.svm.SVC]]:
    return [
        (cutoff, sklearn.base.clone(model))
        for cutoff, model in itertools.product(
            FeatureExtractor.putative_cutoffs(average_cut_position),
            get_putative_svm_models(),
        )
    ]
