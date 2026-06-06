from sklearn.neural_network import MLPClassifier
from sklearn.feature_extraction.text import TfidfVectorizer
from sklearn.model_selection import train_test_split
from tqdm import tqdm
import Bio.SeqUtils.ProtParamData
from Bio.SeqUtils.ProtParam import ProteinAnalysis
import importlib
import itertools
import re
from dataclasses import dataclass
from typing import Iterable, Iterator, Any, Protocol, cast
import sklearn.metrics
import logomaker
import numpy as np
import polars as pl
import sklearn.metrics
import plotly.express as px
import plotly.graph_objects as go
import plotly.io as pio
import contextlib
import io
from IPython.display import clear_output
import src.data_collection; importlib.reload(src.data_collection); clear_output()  # noqa: E703, E702, E402 # fmt: skip
import src.logo_generator; importlib.reload(src.logo_generator); clear_output()  # noqa: E703, E702, E402 # fmt: skip
import src.utils.AdditionalProtParamData; importlib.reload(src.utils.AdditionalProtParamData); clear_output()  # noqa: E703, E702, E402 # fmt: skip
import src.methods.von_heijne; importlib.reload(src.methods.von_heijne);  # noqa: E703, E702, E402 # fmt: skip
import src.processing; importlib.reload(src.processing); clear_output()  # noqa: E703, E702, E402 # fmt: skip
import src.graphics; importlib.reload(src.graphics); clear_output()  # noqa: E703, E702, E402 # fmt: skip
import src.feature_extraction; importlib.reload(src.feature_extraction); clear_output()  # noqa: E703, E702, E402 # fmt: skip
import src.methods.svm; importlib.reload(src.methods.svm); clear_output()  # noqa: E703, E702, E402 # fmt: skip
import src.methods.mlp; importlib.reload(src.methods.mlp); clear_output()  # noqa: E703, E702, E402 # fmt: skip
import src.methods.model_selection; importlib.reload(src.methods.model_selection); clear_output()  # noqa: E703, E702, E402 # fmt: skip
import src.methods.von_heijne; importlib.reload(src.methods.von_heijne); clear_output()  # noqa: E703, E702, E402 # fmt: skip
import src.feature_extraction; importlib.reload(src.feature_extraction); clear_output()  # noqa: E703, E702, E402 # fmt: skip
import src.feature_importance; importlib.reload(src.feature_importance); clear_output()  # noqa: E703, E702, E402 # fmt: skip

# collect data
all_positive, all_negative, positive, negative = src.data_collection.collect_data()
import ankh
import torch

model, tokenizer = ankh.load_base_model()

protein_sequences = [
    "MKALCLLLLPVLGLLVSSKTLCSMEEAINERIQEVAGSLIFRAISSIGLECQSVTSRGDLATCPRGFAVTGCTCGSACGSWDVRAETTCHCQCAGMDWTGARCCRVQPLEHHHHHH",
    "GSHMSLFDFFKNKGSAATATDRLKLILAKERTLNLPYMEEMRKEIIAVIQKYTKSSDIHFKTLDSNQSVETIEVEIILPR",
]

protein_sequences = [list(seq) for seq in protein_sequences]

outputs = tokenizer(
    protein_sequences,
    add_special_tokens=True,
    padding=True,
    is_split_into_words=True,
    return_tensors="pt",
)
with torch.no_grad():
    embeddings = model(
        input_ids=outputs["input_ids"], attention_mask=outputs["attention_mask"]
    )
