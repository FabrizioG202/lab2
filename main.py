import importlib
import itertools
import re
from dataclasses import dataclass
from typing import Iterable, Iterator

import logomaker
import numpy as np
import polars as pl
import sklearn.metrics

import src.data_collection
import src.logo_generator
import src.utils.AdditionalProtParamData

importlib.reload(src.utils.AdditionalProtParamData)
importlib.reload(src.data_collection)

collector = src.data_collection.DataCollector(
    positive_query="""
    (
        (taxonomy_id:2759)
        AND (reviewed:true)
        AND (ft_signal_exp:*)
        AND (fragment:false)
        AND (length:[40 TO *])
        AND (existence:1)
    )
    """,
    negative_query="""
    (
        (reviewed:true)
        AND (fragment:false)
        AND (taxonomy_id:2759)
        AND (length:[40 TO *])
        AND (existence:1)
        AND NOT (ft_signal:*)
        AND (
            (cc_scl_term_exp:SL-0191)
            OR (cc_scl_term_exp:SL-0204)
            OR (cc_scl_term_exp:SL-0039)
            OR (cc_scl_term_exp:SL-0091)
            OR (cc_scl_term_exp:SL-0209)
            OR (cc_scl_term_exp:SL-0173)
        )
    )
    """,
)


collector.clean_cache()

# Setup the working director for the project
# (creates .data, and .imgs)
collector.setup_wd()


# Get the positive examples and cluster them
all_positive = collector.cluster_df(collector.get_positive_examples())

# Get the negative examples and cluster them
all_negative = collector.cluster_df(collector.get_negative_examples())


# positive are only the ones where accession matches cluster_id, so we have one representative per cluster
positive = all_positive.filter(pl.col("accession") == pl.col("cluster_id"))
negative = all_negative.filter(pl.col("accession") == pl.col("cluster_id"))

# Add a column with the sequence neighbouring the Cleavage site.
# ┌───────────┬───────────────────────┬──────────────────────────┐
# │           │                       │                          │
# │           ▼                       ▼                          │
# │LPNTGRLAGCTVFITGASRGIGKAIALKAAKDGANIVIAAKTAQPHPKLLGTIYTAAEEIEA│
# │           ─────────────────────┬───                          │
# │                                ▲                             │
# └────────────────────────────────┼─────────────────────────────┘
#                           cleavage site

K_RESIDUES_BEFORE = 13
K_RESIDUES_AFTER = 2

# K_RESIDUES_BEFORE = 6
# K_RESIDUES_AFTER = 2

# Add a column to the positive which has the motifs (the K residues before and after the cleavage site)
# We assume this does not throw since we selected for residues with SP cleavage site > 14 and sequence length > 90
positive = positive.with_columns(
    pl.concat_str(
        [
            pl.col("sequence").str.slice(
                pl.col("cleavage_site") - K_RESIDUES_BEFORE, K_RESIDUES_BEFORE
            ),
            pl.col("sequence").str.slice(pl.col("cleavage_site"), K_RESIDUES_AFTER),
        ]
    ).alias("motif")
)


# generate the logo
fig = src.logo_generator.generate_logo(
    positive["motif"].to_list(),
    left_bound=K_RESIDUES_BEFORE,
    right_bound=K_RESIDUES_AFTER,
)
fig.savefig(".imgs/positive_logo.svg")

# possible aas
ALPHABET = list("GAVPLIMFWYSTCNQHDEKR")

# split 80% train 20% test
positive_train = positive.sample(fraction=0.8, seed=42)

class PSWM:
    data: np.ndarray
    alphabet: list[str]

    def __init__(self, data : np.ndarray, alphabet: list[str]):
        self.data = data
        self.alphabet = alphabet


    # Compute the pswm for a given list of motifs
    @staticmethod
    def compute_for(motifs : list[str],  alphabet: list[str], motif_length: int | None = None,) -> "PSWM":

        if (motif_length is None):
            motif_length = len(motifs[0])

        # Create PSWM (Position Score Weight Matrix) as a 2D array with shape (len(ALPHABET), K_RESIDUES_BEFORE + K_RESIDUES_AFTER)
        pswm = np.ones((motif_length, len(alphabet)), dtype=float)

        # loop over motifs
        for motif in motifs:
            for i, aa in enumerate(motif):
                if aa in ALPHABET:
                    # TODO: This is really inefficient, we should create a map instead.
                    aa_index = ALPHABET.index(aa)
                    pswm[i, aa_index] += 1


        return PSWM(data= np.log2(
            # Normalize by the number of sequences and 20 pseudocounts
            (pswm / (len(motifs) + 20))
            / [
                src.utils.AdditionalProtParamData.swissprot_composition_2[aa] for aa in ALPHABET
            ],
        ), alphabet=alphabet)



    # Gets the list of scores for all the windows of pswm's size
    # over the sequence.
    def get_scores(self, sequence : str, ):
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





# def cross_validation(folds: int, pwsm: np.ndarray, sequence_data: pl.DataFrame):
sequence_data = pl.concat(
    [
        positive.with_columns(pl.lit(1).alias("label")),  # add a labe with value of 0
        negative.with_columns(pl.lit(0).alias("label")),  # add a labe with value of 1
    ],
    how="diagonal",
)

folds = 5


@dataclass
class _FoldDesc:
    training: list[int]
    testing: list[int]
    validation: list[int]

    @staticmethod
    def generate_all(
        k: int,
        training_count: int = 3,
        testing_count: int = 1,
        validation_count: int = 1,
    ) -> "Iterator[_FoldDesc]":
        assert training_count + testing_count + validation_count == k

        # generate all posible orderings of items in the set {0, 1, ..., k-1}
        for ordering in itertools.permutations(range(k)):
            yield _FoldDesc(
                training=list(ordering[:training_count]),
                testing=list(ordering[training_count : training_count + testing_count]),
                validation=list(ordering[training_count + testing_count :]),
            )

def compute_pswm_threshold(folds: int):




# add a fold column to the dataframe data
sequence_data = sequence_data.with_columns(
    pl.col("accession").map_elements(lambda x: hash(x) % folds).alias("fold")
)

desc = next(_FoldDesc.generate_all(folds))


y_validation = sequence_data.filter(pl.col("fold").is_in(desc.validation))[
    "label"
].to_list()


y_validation_scores = [
    np.max(get_scores_for_sequence(PWSM, s, ALPHABET))
    for s in sequence_data.filter(pl.col("fold").is_in(desc.validation))[
        "sequence"
    ].to_list()
]

precision, recall, thresholds = sklearn.metrics.precision_recall_curve(
    y_validation, y_validation_scores
)

fscore = (
    2
    * (precision * recall)
    / (
        precision + recall + 1e-8  # to avoid division by zero
    )
)
index = np.argmax(fscore)
optimal_threshold = thresholds[index]

y_test_prediction = [
    int(np.max(get_scores_for_sequence(PWSM, s, ALPHABET)) > optimal_threshold)
    for s in sequence_data.filter(pl.col("fold").is_in(desc.testing))[
        "sequence"
    ].to_list()
]

y_test_truth = sequence_data.filter(pl.col("fold").is_in(desc.testing))[
    "label"
].to_list()

confusion_matrix = sklearn.metrics.confusion_matrix(y_test_truth, y_test_prediction)
TN, FP, FN, TP = confusion_matrix.ravel()

accuracy = (TP + TN) / np.sum(confusion_matrix)
sensitivity = TP / (TP + FN) if (TP + FN) > 0 else 0
