# See:
#     https://chanzuckerberg.github.io/cellxgene-census/notebooks/analysis_demo/comp_bio_explore_and_load_lung_data.html
#     http://localhost:8889/notebooks/python_raw/get_dataset.ipynb
#     https://nsforest.readthedocs.io/en/latest/tutorial.html

# TODO: Use pathlib
import logging
from multiprocessing.pool import Pool
import os
import pickle
import re
import requests
import subprocess
from time import sleep
from urllib import parse

from bs4 import BeautifulSoup
import cellxgene_census
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc

import nsforest as ns
from nsforest import nsforesting


DATA_DIR = "data"

CELLXGENE_DOMAIN_NAME = "cellxgene.cziscience.com"
CELLXGENE_API_URL_BASE = f"https://api.{CELLXGENE_DOMAIN_NAME}"
CELLXGENE_DIR = f"{DATA_DIR}/cellxgene"

NSFOREST_DIR = f"{DATA_DIR}/nsforest"

EUTILS_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
EMAIL = "raymond.leclair@gmail.com"
NCBI_API_KEY = os.environ.get("NCBI_API_KEY")
NCBI_API_SLEEP = 1
PUBMED = "pubmed"
PUBMEDCENTRAL = "pmc"

ONTOGPT_DIR = f"{DATA_DIR}/ontogpt"

NCBI_CELL_DIR = f"{DATA_DIR}/ncbi-cell"


def get_lung_obs_and_datasets():

    lung_obs_parquet = f"{NCBI_CELL_DIR}/lung_obs.parquet"
    lung_datasets_parquet = f"{NCBI_CELL_DIR}/lung_datasets.parquet"

    if not os.path.exists(lung_obs_parquet) or not os.path.exists(
        lung_datasets_parquet
    ):

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

        print("Write lung obs parquet")
        lung_obs.to_parquet(lung_obs_parquet)

        print("Finding lung datasets")
        lung_datasets = datasets[datasets["dataset_id"].isin(lung_obs["dataset_id"])]

        print("Write lung datasets parquet")
        lung_datasets.to_parquet(lung_datasets_parquet)

    else:

        print("Read lung obs parquet")
        lung_obs = pd.read_parquet(lung_obs_parquet)

        print("Read lung datasets parquet")
        lung_datasets = pd.read_parquet(lung_datasets_parquet)

    return lung_obs, lung_datasets


def get_title(citation):

    citation_url = None
    title = None

    p1 = re.compile("Publication: (.*) Dataset Version:")
    p2 = re.compile("articleName : '(.*)',")

    selectors = [
        "h1.c-article-title",
        "h1.article-header__title.smaller",
        "div.core-container h1",
        "h1.content-header__title.content-header__title--xx-long",
        "h1#page-title.highwire-cite-title",
    ]

    m1 = p1.search(citation)
    if not m1:
        logging.warning(f"Could not find citation URL for {citation}")
        return citation_url, title

    citation_url = m1.group(1)

    print(f"Getting title for citation URL: {citation_url}")

    response = requests.get(citation_url)

    try_wget = True
    if response.status_code == 200:
        html_data = response.text

        fullsoup = BeautifulSoup(html_data, features="lxml")
        for selector in selectors:
            selected = fullsoup.select(selector)
            if selected:
                if len(selected) > 1:
                    logging.warning(
                        f"Selected more than one element using {selector} on soup from {citation_url}"
                    )
                title = selected[0].text
                try_wget = False
                break

    if try_wget:

        completed_process = subprocess.run(
            ["curl", "-L", citation_url], capture_output=True
        )

        html_data = completed_process.stdout
        fullsoup = BeautifulSoup(html_data, features="lxml")

        found = fullsoup.find_all("script")
        if found and len(found) > 4:
            m2 = p2.search(found[4].text)
            if m2:
                title = m2.group(1)

    print(f"Found title: '{title}' for citation URL: {citation_url}")

    return title


def get_titles(lung_datasets):

    titles = []

    titles_pickle = f"{NCBI_CELL_DIR}/titles.pickle"

    if not os.path.exists(titles_pickle):

        print("Getting titles")
        citations = [c for c in lung_datasets.citation]
        with Pool(8) as p:
            titles = p.map(get_title, citations)
        titles = list(set([title for title in titles if title is not None]))

        print("Dumping titles")
        with open(titles_pickle, "wb") as f:
            pickle.dump(titles, f, pickle.HIGHEST_PROTOCOL)

    else:

        print("Loading titles")
        with open(titles_pickle, "rb") as f:
            titles = pickle.load(f)

    return titles


def get_pmid_for_title(title):

    print(f"Getting PMID for title: '{title}'")

    pmid = None

    search_url = EUTILS_URL + "esearch.fcgi"

    params = {
        "db": PUBMED,
        "term": title,
        "field": "title",
        "retmode": "json",
        # "retmax": 0,
        "email": EMAIL,
        "api_key": NCBI_API_KEY,
    }

    sleep(NCBI_API_SLEEP)

    response = requests.get(search_url, params=parse.urlencode(params, safe=","))

    if response.status_code == 200:
        data = response.json()

        resultcount = int(data["esearchresult"]["count"])
        if resultcount > 1:
            logging.warning(f"PubMed returned more than one result for title: {title}")
            for _pmid in data["esearchresult"]["idlist"]:
                _title = get_title_for_pmid(_pmid)
                if (
                    _title == title + "."
                ):  # PubMedCentral includes period in title, PubMed does not
                    pmid = _pmid
                    print(f"Found PMID: {pmid} for title: '{title}'")
                    break

            if not pmid:
                pmid = data["esearchresult"]["idlist"][0]
                print(f"Using first PMID: {pmid} for title '{title}'")

        else:
            pmid = data["esearchresult"]["idlist"][0]
            print(f"Found PMID: {pmid} for title: '{title}'")

    elif response.status_code == 429:
        logging.error("Too many requests to NCBI API. Try again later, or use API key.")

    else:
        logging.error("Encountered error in searching PubMed: {response.status_code}")

    return pmid


def get_title_for_pmid(pmid):

    title = None

    fetch_url = EUTILS_URL + "efetch.fcgi"

    params = {
        "db": PUBMED,
        "id": pmid,
        "rettype": "xml",
        "email": EMAIL,
        "api_key": NCBI_API_KEY,
    }

    sleep(NCBI_API_SLEEP)

    response = requests.get(fetch_url, params=parse.urlencode(params, safe=","))

    if response.status_code == 200:
        xml_data = response.text

        fullsoup = BeautifulSoup(xml_data, "xml")
        found = fullsoup.find("ArticleTitle")
        if found:
            title = found.text

    else:
        logging.error(
            f"Encountered error in fetching from PubMed: {response.status_code}"
        )

    return title


def get_pmids(titles):

    pmids = []

    pmids_pickle = f"{NCBI_CELL_DIR}/pmids.pickle"

    if not os.path.exists(pmids_pickle):

        print("Getting PMIDs")
        with Pool(8) as p:
            pmids = p.map(get_pmid_for_title, titles)
        pmids = list(set([pmid for pmid in pmids if pmid is not None]))

        print("Dumping PMIDs")
        with open(pmids_pickle, "wb") as f:
            pickle.dump(pmids, f, pickle.HIGHEST_PROTOCOL)

    else:

        print("Loading PMIDs")
        with open(pmids_pickle, "rb") as f:
            pmids = pickle.load(f)

    return pmids


def run_ontogpt_pubmed_annotate(pmid):
    """
    run_ontogpt_pubmed_annotate("38540357")
    """
    output_filename = f"{pmid}.out"
    output_filepath = f"{ONTOGPT_DIR}/{output_filename}"
    if not os.path.exists(output_filepath):
        print(f"Running ontogpt pubmed-annotate for PMID: {pmid}")

        subprocess.run(
            [
                "ontogpt",
                "pubmed-annotate",
                "--template",
                "cell_type",
                pmid,
                "--limit",
                "1",
                "--output",
                output_filepath,
            ],
        )

        print(f"Completed ontogpt pubmed-annotate for PMID: {pmid}")

    else:
        print(f"Ontogpt pubmed-annotate output for PMID: {pmid} exists")


def run_ontogpt(pmids):

    with Pool(8) as p:
        p.map(run_ontogpt_pubmed_annotate, pmids)


def get_dataset(dataset_series):

    collection_id = dataset_series.collection_id
    dataset_id = dataset_series.dataset_id

    dataset_url = f"{CELLXGENE_API_URL_BASE}/curation/v1/collections/{collection_id}/datasets/{dataset_id}"

    response = requests.get(dataset_url)
    response.raise_for_status()

    if response.status_code != 200:
        logging.error(f"Could not get dataset for id {dataset_id}")
        return

    data = response.json()
    if dataset_id != data["dataset_id"]:
        logging.error(
            f"Response dataset id: {data['dataset_id']} does not equal specified dataset id: {dataset_id}"
        )
        return

    assets = data["assets"]
    for asset in assets:

        if asset["filetype"] != "H5AD":
            continue

        dataset_filename = f"{dataset_id}.{asset['filetype']}"
        dataset_filepath = f"{CELLXGENE_DIR}/{dataset_filename}"
        if not os.path.exists(dataset_filepath):

            print(f"Downloading dataset file: {dataset_filepath}")

            with requests.get(asset["url"], stream=True) as response:
                response.raise_for_status()

                with open(dataset_filepath, "wb") as df:
                    for chunk in response.iter_content(chunk_size=1024 * 1024):
                        df.write(chunk)

            print(f"Dataset file: {dataset_filepath} downloaded")

        else:

            print(f"Dataset file: {dataset_filepath} exists")

    return dataset_filename


def get_datasets(datasets_df):

    datasets_series = [row for index, row in datasets_df.iterrows()]
    with Pool(8) as p:
        dataset_filenames = p.map(get_dataset, datasets_series)

    return dataset_filenames


def run_nsforest_on_file(h5ad_filename, cluster_header="cell_type_ontology_term_id"):
    """
    Notes:

    - Some datasets have multiple annotations per sample
    (ex. "broad_cell_type" and "granular_cell_type"). NS-Forest can be
    run on multiple `cluster_header`'s. Combining the parent and child
    markers may improve classification results.

    - `adata.var_names` must be unique. If there is a problem, usually
    it can be solved by assigning `adata.var.index =
    adata.var["ensembl_id"]`.

    - Some datasets are too large and need to be downsampled to be run
    through the pipeline. When downsampling, be sure to have all the
    granular cluster annotations represented.

    - Only run ns.pp.dendrogram() if there is no pre-defined dendrogram
    order. This step can still be run with no effects, but the runtime
    may increase.
    """
    # Assign results directory
    results_dirname = h5ad_filename.split(".")[0]
    results_dirpath = f"{NSFOREST_DIR}/{results_dirname}"
    if not os.path.exists(results_dirpath):
        os.mkdir(results_dirpath)

        print(f"Loading unprocessed AnnData file: {h5ad_filename}")
        h5ad_filepath = f"{CELLXGENE_DIR}/{h5ad_filename}"
        up_adata = sc.read_h5ad(h5ad_filepath)

        print("Generating scanpy dendrogram")
        # Dendrogram order is stored in
        # `pp_adata.uns["dendrogram_cluster"]["categories_ordered"]`
        pp_adata = up_adata.copy()
        pp_adata.obs[cluster_header] = pp_adata.obs[cluster_header].astype(str)
        pp_adata.obs[cluster_header] = pp_adata.obs[cluster_header].astype("category")
        pp_adata = ns.pp.dendrogram(
            pp_adata,
            cluster_header,
            save=False,
            output_folder=results_dirpath,
            outputfilename_suffix=cluster_header,
        )

        print("Calculating cluster medians per gene")
        pp_adata = ns.pp.prep_medians(pp_adata, cluster_header)

        print("Calculating binary scores per gene per cluster")
        pp_adata = ns.pp.prep_binary_scores(pp_adata, cluster_header)

        pp_h5ad_filename = f"pp_{h5ad_filename}"
        pp_h5ad_filepath = f"{results_dirpath}/{pp_h5ad_filename}"
        print(f"Saving preprocessed AnnData file: {pp_h5ad_filepath}")
        pp_adata.write_h5ad(pp_h5ad_filepath)

        print(f"Running NS-Forest for preprocessed AnnData file: {pp_h5ad_filename}")
        results = nsforesting.NSForest(
            pp_adata,
            cluster_header,
            output_folder=f"{results_dirpath}/",
            outputfilename_prefix=cluster_header,
        )

    else:
        print(f"Completed NS-Forest for unprocessed AnnData file: {h5ad_filename}")


def main():

    lung_obs, lung_datasets = get_lung_obs_and_datasets()
    titles = get_titles(lung_datasets)
    pmids = get_pmids(titles)
    run_ontogpt(pmids)
    dataset_filenames = get_datasets(lung_datasets)
    for dataset_filename in dataset_filenames:
        try:
            run_nsforest_on_file(dataset_filename)
        except Exception as ex:
            print(f"Could not run NS-Forest for unprocessed AnnData file: {dataset_filename}")


if __name__ == "__main__":
    __spec__ = None  # Workaround for Pool() in IPython
    main()
