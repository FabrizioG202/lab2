import logging
import os
import re
import subprocess
import tempfile
import warnings
from typing import Self

import polars as pl
import requests
from Bio import SeqIO
from Bio.Align import parse
from debugpy.launcher.debuggee import process
from numpy import positive
from parso.python.tree import Literal
from requests.adapters import HTTPAdapter, Retry
from sklearn.model_selection import train_test_split
from tqdm import tqdm

from src.run_mmseqs import run_mmseqs


def query_uniprot_streamed(query: str, fields: list[str], file_name: str):
    # Code below is directly from uniprot's own API docs.
    re_next_link = re.compile(r'<(.+)>; rel="next"')
    retries = Retry(total=5, backoff_factor=0.25, status_forcelist=[500, 502, 503, 504])
    session = requests.Session()
    session.mount("https://", HTTPAdapter(max_retries=retries))

    def __get_next_link(headers):
        if "Link" in headers:
            match = re_next_link.match(headers["Link"])
            if match:
                return match.group(1)

    def __get_batch(batch_url):
        while batch_url:
            response = session.get(batch_url)
            response.raise_for_status()
            total = response.headers["x-total-results"]
            yield response, total
            batch_url = __get_next_link(response.headers)

    fields_data = ",".join(fields)

    query_encoded = requests.utils.quote(query)

    url = f"https://rest.uniprot.org/uniprotkb/search?format=tsv&query={query_encoded}&size=500&fields={fields_data}"

    with open(file_name, "w", encoding="utf8") as file:
        total = None

        # header line
        file.write("\t".join(fields) + "\n")

        with tqdm(total=total, desc="Fetching UniProt data") as pbar:
            for batch, total in __get_batch(url):
                # Discard the header line, we use the fields
                lines = batch.text.splitlines()[1:]

                file.write("\n".join(lines) + "\n")
                pbar.update(len(lines))
                pbar.total = int(total)


def query_uniprot_json_objs(query: str, page_size: int = 500):
    # Code below is directly from uniprot's own API docs.

    re_next_link = re.compile(r'<(.+)>; rel="next"')
    retries = Retry(total=5, backoff_factor=0.25, status_forcelist=[500, 502, 503, 504])
    session = requests.Session()
    session.mount("https://", HTTPAdapter(max_retries=retries))

    def __get_next_link(headers):
        if "Link" in headers:
            match = re_next_link.match(headers["Link"])
            if match:
                return match.group(1)

    def __get_batch(batch_url):
        while batch_url:
            response = session.get(batch_url)
            response.raise_for_status()
            total = response.headers["x-total-results"]
            yield response, total
            batch_url = __get_next_link(response.headers)

    query_encoded = requests.utils.quote(query)

    url = f"https://rest.uniprot.org/uniprotkb/search?format=json&query={query_encoded}&size={page_size}"

    with tqdm(desc="Fetching UniProt data") as pbar:
        for batch, total in __get_batch(url):
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
            logging.info(
                "Using cached positive examples. Set ignore_cache=True to fetch new data from UniProt."
            )
            return pl.read_csv(cache_file, separator="\t")

        # Fetch the data from accession
        # TODO: This can be made a bit better by using a function and a generator I think
        rows = []

        for o in query_uniprot_json_objs(self.positive_query):
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
            logging.info(
                "Using cached negative examples. Set ignore_cache=True to fetch new data from UniProt."
            )
            return pl.read_csv(cache_file, separator="\t")

        rows = []

        # Fetch the data from accession
        # TODO: This can be made a bit better by using a function and a generator I think
        for o in query_uniprot_json_objs(self.negative_query):
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

    def cluster_df(self: Self, df: pl.DataFrame) -> pl.DataFrame:
        """
        Cluster the sequences in the dataframe using mmseqs2 and add a column with the cluster id.
        """
        with tempfile.NamedTemporaryFile(
            "w",
            delete=True,
            encoding="utf8",
            dir=".data",
        ) as tmp_fasta:
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
            # TODO: throw if these do not exist...
            os.remove(tmp_fasta.name + "_clustered_cluster.tsv")
            os.remove(tmp_fasta.name + "_clustered_all_seqs.fasta")
            os.remove(tmp_fasta.name + "_clustered_rep_seq.fasta")

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
