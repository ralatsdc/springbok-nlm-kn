# https://chanzuckerberg.github.io/cellxgene-census/notebooks/analysis_demo/comp_bio_explore_and_load_lung_data.html
import logging
from multiprocessing.pool import ThreadPool
import os
import subprocess

import cellxgene_census
import pandas as pd


def get_lung_datasets():

    lung_parquet = "lung.parquet"
    if not os.path.exists(lung_parquet):

        print("Opening soma")
        census = cellxgene_census.open_soma(census_version="latest")

        print("Collecting all datasets")
        datasets = census["census_info"]["datasets"].read().concat().to_pandas()

        print("Collecting lung obs")
        lung_obs = (
            census["census_data"]["homo_sapiens"]
            .obs.read(
                value_filter="tissue_general == 'lung' and is_primary_data == True"
            )
            .concat()
            .to_pandas()
        )

        census.close()

        print("Finding lung datasets")
        lung_datasets = datasets[datasets["dataset_id"].isin(lung_obs["dataset_id"])]

        print("Write lung datasets parquet")
        lung_datasets.to_parquet(lung_parquet)

    else:

        print("Read lung datasets parquet")
        lung_datasets = pd.read_parquet(lung_parquet)

    return lung_datasets


def run_ontogpt_pubmed_annotate(pmid):
    subprocess.run(
        [
            "ontogpt",
            "pubmed-annotate",
            "--template",
            "cell_type",
            f"{pmid}",
            "--limit",
            "1",
            "--output",
            f"{pmid}.out",
        ],
        capture_output=True,
    )


if __name__ == "__main__":
    run_ontogpt_pubmed_annotate("38540357")

    """
    with ThreadPool(3) as p:
        result = p.map(f, [1, 2, 3, 4, 5, 6])
    """
