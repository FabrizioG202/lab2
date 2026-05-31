import contextlib
import io
from typing import Any

import Bio.SeqUtils.ProtParamData
import numpy as np
from Bio.SeqUtils.ProtParam import ProteinAnalysis

import src.utils.AdditionalProtParamData


# A function to get the compositions in aa for the given sequence
# TODO: Make the alphabet a parameter...
def get_composition(sequence: str, alphabet: list[str]) -> dict[str, float]:
    composition = {aa: 0 for aa in alphabet}
    for aa in sequence:
        if aa in composition:
            composition[aa] += 1
    total = len(sequence)
    return {aa: count / total for aa, count in composition.items()}


def get_all_features(sequence: str, alphabet: list[str]) -> dict[str, float]:
    window_size = 5

    # pad the sequence with X so that with the given
    # window size, we can always get a full window of residues around the cleavage site.
    def _pad_sequence(seq: str) -> str:
        padding = "X" * (window_size // 2)
        return padding + seq + padding

    def _calc_params(name: str, values: list[float]) -> dict[str, int | float | Any]:
        return {
            f"max_{name}": np.max(values),
            f"avg_{name}": np.mean(values),
        }

    with contextlib.redirect_stderr(io.StringIO()):
        analysis = ProteinAnalysis(_pad_sequence(sequence))

        out = {
            **{
                f"comp_{aa}": v
                for aa, v in get_composition(sequence[:22], alphabet).items()
            },
            **_calc_params(
                "hydrophobicity",
                analysis.protein_scale(Bio.SeqUtils.ProtParamData.kd, 5),
            ),
            **_calc_params(
                "alpha_helix_tendency",
                analysis.protein_scale(
                    src.utils.AdditionalProtParamData.alpha_helix_tendency, 5
                ),
            ),
            **_calc_params(
                "transmembrane_tendency",
                analysis.protein_scale(
                    src.utils.AdditionalProtParamData.transmemberane_tendency, 5
                ),
            ),
            **_calc_params(
                "bulkiness",
                analysis.protein_scale(src.utils.AdditionalProtParamData.bulkiness, 5),
            ),
            **_calc_params(
                "charge",
                analysis.protein_scale(src.utils.AdditionalProtParamData.polarity, 5),
            ),
        }

    return out
