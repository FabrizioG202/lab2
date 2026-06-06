import itertools
from dataclasses import dataclass
from typing import Iterator, Union

import numpy as np
import polars as pl
import sklearn.metrics
import tqdm

import src.graphics
import src.utils.AdditionalProtParamData


@dataclass
class _FoldValidationMetric:
    optimal_threshold: float
    TP: int
    TN: int
    FP: int
    FN: int

    @property
    def accuracy(self) -> float:
        return (self.TP + self.TN) / (self.TP + self.TN + self.FP + self.FN)

    @property
    def sensitivity(self) -> float:
        return self.TP / (self.TP + self.FN) if (self.TP + self.FN) > 0 else 0


class PSWM:
    data: np.ndarray
    alphabet: list[str]

    def __init__(self, data: np.ndarray, alphabet: list[str]):
        self.data = data
        self.alphabet = alphabet

    # Compute the pswm for a given list of motifs
    @staticmethod
    def compute_for(
        motifs: list[str],
        alphabet: list[str],
        motif_length: Union[int, None] = None,
    ) -> "PSWM":

        if motif_length is None:
            motif_length = len(motifs[0])

        # Create PSWM (Position Score Weight Matrix) as a 2D array with shape (len(ALPHABET), K_RESIDUES_BEFORE + K_RESIDUES_AFTER)
        pswm = np.ones((motif_length, len(alphabet)), dtype=float)

        # loop over motifs
        for motif in motifs:
            for i, aa in enumerate(motif):
                if aa in alphabet:
                    # TODO: This is really inefficient, we should create a map instead.
                    aa_index = alphabet.index(aa)
                    pswm[i, aa_index] += 1

        return PSWM(
            data=np.log2(
                # Normalize by the number of sequences and 20 pseudocounts
                (pswm / (len(motifs) + 20))
                / [
                    src.utils.AdditionalProtParamData.swissprot_composition[aa]
                    for aa in alphabet
                ],
            ),
            alphabet=alphabet,
        )

    # Gets the list of scores for all the windows of pswm's size
    # over the sequence.
    def get_scores(self, sequence: str):
        scores = []

        window_size = self.data.shape[0]

        # We select only the first 90 aa
        seq = sequence[:90]

        for i in range(len(seq) - window_size + 1):
            window = seq[i : i + window_size]
            score = 0

            for j, aa in enumerate(window):
                if aa in self.alphabet:
                    aa_index = self.alphabet.index(aa)
                    score += self.data[j, aa_index]

            scores.append(score)

        return scores

    def compute_optimal_threshold(
        self, *, folds: int, all_sequences: pl.DataFrame
    ) -> tuple[float, list[_FoldValidationMetric]]:
        assert "accession" in all_sequences.columns
        assert "sequence" in all_sequences.columns
        assert "label" in all_sequences.columns

        history = []

        # add a fold column to the dataframe data
        sequence_data = all_sequences.with_columns(
            pl.col("accession").map_elements(lambda x: hash(x) % folds).alias("fold")
        )

        for desc in tqdm.tqdm(
            _FoldDesc.generate_all(folds),
            desc="Computing optimal threshold for each fold",
        ):
            validation_df = sequence_data.filter(pl.col("fold").is_in(desc.validation))
            test_df = sequence_data.filter(pl.col("fold").is_in(desc.testing))

            # compute the preicsion, recall and thresholds based on the validation set for this fold.
            precision, recall, thresholds = sklearn.metrics.precision_recall_curve(
                validation_df["label"].to_list(),
                [
                    np.max(self.get_scores(seq))
                    for seq in validation_df["sequence"].to_list()
                ],
            )

            fscore = (
                2
                * (precision * recall)
                / (
                    precision + recall + 1e-8  # to avoid division by zero
                )
            )
            index = np.argmax(fscore)

            # this is the optimal threshold
            optimal_threshold = thresholds[index]

            confusion_matrix = sklearn.metrics.confusion_matrix(
                test_df["label"].to_list(),
                [
                    int(np.max(self.get_scores(seq)) >= optimal_threshold)
                    for seq in test_df["sequence"].to_list()
                ],
            )

            TN, FP, FN, TP = confusion_matrix.ravel()

            history.append(
                _FoldValidationMetric(
                    optimal_threshold=optimal_threshold,
                    TP=TP,
                    TN=TN,
                    FP=FP,
                    FN=FN,
                )
            )

        return np.mean([h.optimal_threshold for h in history]), history


@dataclass
class _FoldDesc:
    training: set[int]
    testing: set[int]
    validation: set[int]

    @staticmethod
    def generate_all(
        k: int,
        training_count: int = 3,
        testing_count: int = 1,
        validation_count: int = 1,
    ) -> "set[_FoldDesc]":
        assert training_count + testing_count + validation_count == k

        # This is horribly inefficient, since we are converting the stuff
        # back to a set after we computed all, but for small fold values
        # this should suffice
        out = []

        # generate all posible orderings of items in the set {0, 1, ..., k-1}
        for ordering in itertools.permutations(range(k)):
            out.append(
                _FoldDesc(
                    training=set(ordering[:training_count]),
                    testing=set(
                        ordering[training_count : training_count + testing_count]
                    ),
                    validation=set(ordering[training_count + testing_count :]),
                )
            )

        return set(out)

    def __hash__(self):
        return hash(
            (
                frozenset(self.training),
                frozenset(self.testing),
                frozenset(self.validation),
            )
        )

    def __eq__(self, other):
        if not isinstance(other, _FoldDesc):
            return False

        return (
            frozenset(self.training) == frozenset(other.training)
            and frozenset(self.testing) == frozenset(other.testing)
            and frozenset(self.validation) == frozenset(other.validation)
        )


def compute_pswm(
    combined_train: pl.DataFrame, combined_test: pl.DataFrame
) -> tuple[PSWM, src.graphics.ConfusionMatrix]:
    # Compute the PSWM on the training set, using only the positive examples
    # with a motif (i.e. where we could extract the motif).
    pswm = PSWM.compute_for(
        [m for m in combined_train["motif"].to_list() if m is not None],
        alphabet=list("GAVPLIMFWYSTCNQHDEKR"),
    )

    # Use CV to compute the optimal threshold for the PSWM scores on the training set.
    threshold, history = pswm.compute_optimal_threshold(
        folds=5,
        all_sequences=combined_train.select(
            [
                "accession",
                "sequence",
                "label",
            ]
        ),
    )

    # Now, evaluate on the test set and compute the confusion matrix for the PSWM method and save it as an image
    return pswm, src.graphics.ConfusionMatrix(
        sklearn.metrics.confusion_matrix(
            combined_test["label"].to_list(),
            [
                1 if max(pswm.get_scores(motif[:90])) >= threshold else 0
                for motif in combined_test["sequence"].to_list()
            ],
        )
    )
