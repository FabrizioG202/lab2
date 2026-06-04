import logging
import os
import re
import subprocess
import tempfile
import warnings
from typing import Callable, Self

import plotly.graph_objects as go
import polars as pl
import requests
from Bio import SeqIO
from Bio.Align import parse
from debugpy.launcher.debuggee import process
from parso.python.tree import Literal
from requests.adapters import HTTPAdapter, Retry
from sklearn.model_selection import train_test_split
from tqdm import tqdm

from src.run_mmseqs import run_mmseqs


def _query_uniprot_json_objs(query: str, page_size: int = 500):
    # Code below is directly from uniprot's own API docs.

    re_next_link = re.compile(r'<(.+)>; rel="next"')
    retries = Retry(total=5, backoff_factor=0.25, status_forcelist=[500, 502, 503, 504])
    session = requests.Session()
    session.mount("https://", HTTPAdapter(max_retries=retries))

    def _get_next_link(headers):
        if "Link" in headers:
            match = re_next_link.match(headers["Link"])
            if match:
                return match.group(1)

    def _get_batch(batch_url):
        while batch_url:
            response = session.get(batch_url)
            response.raise_for_status()
            total = response.headers["x-total-results"]
            yield response, total
            batch_url = _get_next_link(response.headers)

    query_encoded = requests.utils.quote(query)

    url = f"https://rest.uniprot.org/uniprotkb/search?format=json&query={query_encoded}&size={page_size}"

    with tqdm(desc="Fetching UniProt data") as pbar:
        for batch, total in _get_batch(url):
            results = batch.json()["results"]

            pbar.update(len(results))
            pbar.total = int(total)

            yield from results


class DataCollector:
    positive_query: str
    negative_query: str

    def __init__(
        self,
        positive_query: str,
        negative_query: str,
    ) -> None:
        self.positive_query = positive_query
        self.negative_query = negative_query

    def get_positive_examples(
        self: Self, ignore_cache: bool = False, cache_file: str = ".data/positive.tsv"
    ) -> pl.DataFrame:
        """
        Get the positive examples from UniProt.
        If ignore_cache is False, it will check if the data is already cached in `cache_file`
        and use that instead of fetching new data from UniProt.
        """

        # Check if cached data exists
        if not ignore_cache and os.path.exists(cache_file):
            print(
                "Using cached positive examples. Set ignore_cache=True to fetch new data from UniProt."
            )
            return pl.read_csv(cache_file, separator="\t")

        # Fetch the data from accession
        # TODO: This can be made a bit better by using a function and a generator I think
        rows = []

        objects = list(_query_uniprot_json_objs(self.positive_query))

        print(f"Total Number of positive examples fetched from UniProt: {len(objects)}")

        for o in objects:
            features = o["features"]
            signal_features = [f for f in features if f["type"] == "Signal"]

            if len(signal_features) == 0:
                continue

            signal_end = signal_features[0]["location"]["end"]["value"]

            if not signal_end:
                continue

            signal_start = signal_features[0]["location"]["start"]["value"]

            # filter out SP shorter than 14 residues
            if signal_end - signal_start < 14:
                continue

            rows.append(
                {
                    "accession": o["primaryAccession"],
                    "organism_name": o["organism"]["scientificName"],
                    "cleavage_site": signal_end,
                    "kingdom": o["organism"]["lineage"][1],
                    "protein_length": o["sequence"]["length"],
                    "sequence": o["sequence"]["value"],
                }
            )

        # Cache the data
        df = pl.DataFrame(rows)
        df.write_csv(cache_file, separator="\t")

        return df

    def get_negative_examples(
        self: Self, ignore_cache: bool = False, cache_file: str = ".data/negative.tsv"
    ) -> pl.DataFrame:
        """
        Get the negative examples from UniProt.
        If ignore_cache is False, it will check if the data is already cached in `cache_file`
        and use that instead of fetching new data from UniProt.
        """

        # Check if cached data exists
        if not ignore_cache and os.path.exists(cache_file):
            print(
                "Using cached negative examples. Set ignore_cache=True to fetch new data from UniProt."
            )
            return pl.read_csv(cache_file, separator="\t")

        rows = []

        # Fetch the data from accession
        # TODO: This can be made a bit better by using a function and a generator I think
        for o in _query_uniprot_json_objs(self.negative_query):
            # TODO: I think this can be inlined
            # Wether the protein has a transmembrane domain in the first 90 residues
            has_transmembrane = False
            for feature in o["features"]:
                if (
                    feature["type"] == "Transmembrane"
                    and "Helical" in feature["description"]
                    and feature["location"]["start"]["value"] <= 90
                ):
                    has_transmembrane = True

            rows.append(
                {
                    "accession": o["primaryAccession"],
                    "organism_name": o["organism"]["scientificName"],
                    "has_transmembrane": has_transmembrane,
                    "kingdom": o["organism"]["lineage"][1],
                    "protein_length": o["sequence"]["length"],
                    "sequence": o["sequence"]["value"],
                }
            )

        # Cache the data
        df = pl.DataFrame(rows)
        df.write_csv(cache_file, separator="\t")

        return df

    def cluster_df(
        self: Self,
        get_df: Callable[[], pl.DataFrame],
        cache_prefix: str,
        ignore_cache: bool = False,
    ) -> pl.DataFrame:
        """
        Cluster the sequences in the dataframe using mmseqs2 and add a column with the cluster id.
        """
        cache_file = f".data/{cache_prefix}.clustered.tsv"

        # Check if cached data exists
        if not ignore_cache and os.path.exists(cache_file):
            print(
                "Using cached clustered data. Set ignore_cache=True to fetch new data from UniProt."
            )
            return pl.read_csv(cache_file, separator="\t")

        # get the dataframe
        df = get_df()

        tmp_fasta_path = ".data/tmp.fasta"
        with open(tmp_fasta_path, "w", encoding="utf8") as tmp_fasta:
            tmp_fasta.write(
                "\n".join(
                    [
                        f">{r[0]}\n{r[1]}"
                        for r in df.select(["accession", "sequence"]).iter_rows()
                    ]
                )
            )

        # Run mmseqs2 to cluster the sequences at 30% identity
        run_mmseqs(
            tmp_fasta.name,
            tmp_fasta.name + "_clustered",
        )

        # check that the necessary files were generated
        for suffix in [
            "_clustered_cluster.tsv",
            "_clustered_all_seqs.fasta",
            "_clustered_rep_seq.fasta",
        ]:
            if not os.path.exists(tmp_fasta.name + suffix):
                raise FileNotFoundError(
                    f"Expected file {tmp_fasta.name + suffix} was not generated by mmseqs2."
                )

        # get the index file
        index_file = tmp_fasta.name + "_clustered_cluster.tsv"

        # read the index file and create a mapping from accession to cluster id
        cluster_map = {}
        with open(index_file, "r", encoding="utf8") as f:
            for line in f:
                cluster_id, accession = line.strip().split("\t")
                cluster_map[accession.strip()] = cluster_id.strip()

        # add a column to the dataframe with the cluster id
        df = df.with_columns(
            pl.col("accession").replace_strict(cluster_map).alias("cluster_id")
        )

        # remove the temp files
        os.remove(tmp_fasta.name + "_clustered_cluster.tsv")
        os.remove(tmp_fasta.name + "_clustered_all_seqs.fasta")
        os.remove(tmp_fasta.name + "_clustered_rep_seq.fasta")
        os.remove(tmp_fasta.name)

        # save the clustered data to cache
        df.write_csv(cache_file, separator="\t")
        print(f"Clustered data saved to {cache_file}")

        return df

    def setup_wd(self: Self) -> None:
        """
        Setup the working directory for the project by creating the .data and .imgs directories if they do not exist.
        """
        if not os.path.exists(".data"):
            os.makedirs(".data")

        if not os.path.exists(".imgs"):
            os.makedirs(".imgs")

    def clean_cache(self: Self) -> None:
        """
        Clean the data cache directory by removing all files in it.
        """
        for file in os.listdir(".data"):
            file_path = os.path.join(".data", file)
            if os.path.isfile(file_path):
                os.remove(file_path)

    def save_dataset(
        dataframe: pl.DataFrame,
        name: str,
    ) -> None:

        # Save TSV
        dataframe.write_csv(f".data/{name}.tsv", separator="\t")

        # Save FASTA
        with open(f".data/{name}.fa", "w", encoding="utf8") as f:
            f.write(
                "\n".join(
                    [
                        f">{r[0]}\n{r[1]}"
                        for r in dataframe.select(["accession", "sequence"]).iter_rows()
                    ]
                )
            )


def plot_lengths(
    positive: pl.DataFrame,
    negative: pl.DataFrame,
) -> go.Figure:
    positive_lengths = positive["sequence"].str.len_chars().to_list()
    negative_lengths = negative["sequence"].str.len_chars().to_list()

    fig = go.Figure()
    fig.add_trace(
        go.Histogram(
            x=positive_lengths,
            name="Positive",
            marker_color="royalblue",
            opacity=0.65,
            nbinsx=40,
        )
    )
    fig.add_trace(
        go.Histogram(
            x=negative_lengths,
            name="Negative",
            marker_color="indianred",
            opacity=0.65,
            nbinsx=40,
        )
    )

    fig.update_layout(
        title="Sequence Length Distribution",
        xaxis_title="Sequence length",
        yaxis_title="Count",
        barmode="overlay",
    )

    return fig


def plot_kingdom_distribution(
    positive: pl.DataFrame,
    negative: pl.DataFrame,
) -> go.Figure:
    kingdom_counts = (
        pl.concat(
            [
                positive.select(pl.col("kingdom")).with_columns(
                    pl.lit("Positive").alias("dataset")
                ),
                negative.select(pl.col("kingdom")).with_columns(
                    pl.lit("Negative").alias("dataset")
                ),
            ],
            how="diagonal",
        )
        .group_by(["dataset", "kingdom"])
        .len()
        .sort(["dataset", "len"], descending=[False, True])
    )

    fig = go.Figure()
    for dataset_name, color in [("Positive", "royalblue"), ("Negative", "indianred")]:
        subset = kingdom_counts.filter(pl.col("dataset") == dataset_name)
        fig.add_trace(
            go.Pie(
                labels=subset["kingdom"].to_list(),
                values=subset["len"].to_list(),
                name=dataset_name,
                marker=dict(colors=None if dataset_name == "Positive" else None),
                domain=dict(column=0 if dataset_name == "Positive" else 1),
            )
        )

    fig.update_layout(
        title="Kingdom distribution of positive and negative datasets",
        grid=dict(rows=1, columns=2),
        annotations=[
            dict(text="Positive", x=0.20, y=1.08, showarrow=False),
            dict(text="Negative", x=0.80, y=1.08, showarrow=False),
        ],
    )

    return fig


def collect_data(
    purge_cache: bool = False,
) -> tuple[
    pl.DataFrame,
    pl.DataFrame,
    pl.DataFrame,
    pl.DataFrame,
]:
    """Setups the data for the pipeline"""
    collector = DataCollector(
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

    # (if necessary, clean the cache)
    if purge_cache:
        collector.clean_cache()

    # Setup the working director for the project
    # (creates .data, and .imgs)
    collector.setup_wd()

    # These contain information about the clusters, but there still are multiple sequences per cluster.
    all_positive = collector.cluster_df(
        lambda: collector.get_positive_examples(),
        cache_prefix="positive",
    )
    all_negative = collector.cluster_df(
        lambda: collector.get_negative_examples(),
        cache_prefix="negative",
    )

    # positive are only the ones where accession matches cluster_id, so we have one representative per cluster
    positive = all_positive.filter(pl.col("accession") == pl.col("cluster_id"))
    negative = all_negative.filter(pl.col("accession") == pl.col("cluster_id"))

    return all_positive, all_negative, positive, negative
