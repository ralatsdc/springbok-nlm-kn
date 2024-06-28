import logging
from multiprocessing.pool import Pool
import os  # TODO: Use pathlib?
import re
import requests
import subprocess
from time import sleep
from urllib import parse

from bs4 import BeautifulSoup
import cellxgene_census
import pandas as pd
import scanpy as sc

import nsforest as ns
from nsforest import nsforesting


DATA_DIR = "data"

CELLXGENE_DOMAIN_NAME = "cellxgene.cziscience.com"
CELLXGENE_API_URL_BASE = f"https://api.{CELLXGENE_DOMAIN_NAME}"
CELLXGENE_DIR = f"{DATA_DIR}/cellxgene"

HTTPS_SLEEP = 1

NSFOREST_DIR = f"{DATA_DIR}/nsforest-2024-06-27"
TOTAL_COUNTS = 5000  # TODO: Select a more sensible value

EUTILS_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
NCBI_EMAIL = os.environ.get("NCBI_EMAIL")
NCBI_API_KEY = os.environ.get("NCBI_API_KEY")
NCBI_API_SLEEP = 1
PUBMED = "pubmed"
PUBMEDCENTRAL = "pmc"

ONTOGPT_DIR = f"{DATA_DIR}/ontogpt"

NCBI_CELL_DIR = f"{DATA_DIR}/ncbi-cell"


def get_lung_obs_and_datasets():
    """Use the CELLXGENE Census to obtain all unprocessed human lung
    cell metadata and datasets, then write the resulting Pandas
    DataFrames to parquet files, or, if the files exist, read them.

    Parameters
    ----------
    None

    Returns
    -------
    lung_obs : pd.DataFrame
        DataFrame containing unprocessed dataset metadata
    lung_datasets : pd.DataFrame
        DataFrame containing unprocessed dataset descriptions
    """
    # Create and write, or read DataFrames
    lung_obs_parquet = f"{NCBI_CELL_DIR}/up_lung_obs.parquet"
    lung_datasets_parquet = f"{NCBI_CELL_DIR}/up_lung_datasets.parquet"
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

        print("Closing soma")
        census.close()

        print("Writing unprocessed lung obs parquet")
        lung_obs.to_parquet(lung_obs_parquet)

        print("Finding unprocessed lung datasets")
        lung_datasets = datasets[datasets["dataset_id"].isin(lung_obs["dataset_id"])]

        print("Writing unprocessed lung datasets parquet")
        lung_datasets.to_parquet(lung_datasets_parquet)

    else:

        print("Reading unprocessed lung obs parquet")
        lung_obs = pd.read_parquet(lung_obs_parquet)

        print("Reading unprocessed lung datasets parquet")
        lung_datasets = pd.read_parquet(lung_datasets_parquet)

    return lung_obs, lung_datasets


def append_titles_pmids_and_dataset_h5ad_files(up_lung_datasets):
    """Append titles and PubMed identifiers (PMIDs) corresponding to
    the dataset, and the dataset filename to the CELLXGENE dataset
    DataFrame, download the dataset file, and write the preprocessed
    DataFrame to a parquet file, or, if the file exists, read it.

    Parameters
    ----------
    up_lung_datasets : pd.DataFrame
        DataFrame containing unprocessed dataset descriptions

    Returns
    -------
    pp_lung_datasets : pd.DataFrame
        DataFrame containing preprocessed dataset descriptions
    """
    pp_lung_datasets_parquet = f"{NCBI_CELL_DIR}/pp_lung_datasets.parquet"
    if not os.path.exists(pp_lung_datasets_parquet):

        pp_lung_datasets = up_lung_datasets.copy()
        pp_lung_datasets = append_titles(pp_lung_datasets)
        pp_lung_datasets = append_pmids(pp_lung_datasets)
        pp_lung_datasets = append_and_download_dataset_h5ad_files(pp_lung_datasets)

        print("Writing preprocessed lung datasets parquet")
        pp_lung_datasets.to_parquet(pp_lung_datasets_parquet)

    else:

        print("Reading preprocessed lung datasets parquet")
        pp_lung_datasets = pd.read_parquet(pp_lung_datasets_parquet)

    return pp_lung_datasets


def append_titles(lung_datasets):
    """Get and append titles for each dataset citation using a
    subprocess pool.

    Parameters
    ----------
    lung_datasets : pd.DataFrame
        DataFrame containing dataset descriptions

    Returns
    ----------
    lung_datasets : pd.DataFrame
        DataFrame containing dataset descriptions with titles appended
    """
    print("Getting titles")
    with Pool(8) as p:
        titles = p.map(get_title, lung_datasets["citation"])
    lung_datasets["citation_title"] = titles

    return lung_datasets


def append_pmids(lung_datasets):
    """Get and append PMIDs for each citation title using a subprocess
    pool.

    Parameters
    ----------
    lung_datasets : pd.DataFrame
        DataFrame containing dataset descriptions with titles appended

    Returns
    ----------
    lung_datasets : pd.DataFrame
        DataFrame containing dataset descriptions with titles and
        PMIDs appended
    """
    print("Getting PMIDs")
    with Pool(8) as p:
        pmids = p.map(get_pmid_for_title, lung_datasets["citation_title"])
    lung_datasets["citation_pmid"] = pmids

    return lung_datasets


def append_and_download_dataset_h5ad_files(lung_datasets):
    """Get and append dataset filenames for each dataset using a
    subprocess pool.

    Parameters
    ----------
    lung_datasets : pd.DataFrame
        DataFrame containing dataset descriptions

    Returns
    ----------
    lung_datasets : pd.DataFrame
        DataFrame containing dataset descriptions with dataset
        filenames appended
    """
    datasets_series = [row for index, row in lung_datasets.iterrows()]
    print("Getting dataset files")
    with Pool(8) as p:
        dataset_h5ad_files = p.map(get_and_download_dataset_h5ad_file, datasets_series)
    lung_datasets["dataset_h5ad_file"] = dataset_h5ad_files

    return lung_datasets


def get_title(citation):
    """Get the title given a dataset citation. Note that only wget
    succeeded for Cell Press journals, and neither requests nor wget
    succeeded for The EMBO Journal and Science.

    Parameters
    ----------
    citation : str
        Dataset citation

    Returns
    -------
    title : str
        Title of publication associated with the dataset
    """
    # Need a default return value
    title = None

    # Compile patterns for finding the publication URL and article
    # title
    p1 = re.compile("Publication: (.*) Dataset Version:")
    p2 = re.compile("articleName : '(.*)',")

    # Assign CSS selectors for selecting article title elements
    selectors = [
        "h1.c-article-title",
        "h1.article-header__title.smaller",
        "div.core-container h1",
        "h1.content-header__title.content-header__title--xx-long",
        "h1#page-title.highwire-cite-title",
    ]

    # Find the publication URL
    m1 = p1.search(citation)
    if not m1:
        logging.warning(f"Could not find citation URL for {citation}")
        return
    citation_url = m1.group(1)
    print(f"Getting title for citation URL: {citation_url}")

    # Attempt to get the publication page using requests
    print(f"Trying requests")
    sleep(HTTPS_SLEEP)
    response = requests.get(citation_url)
    try_wget = True
    if response.status_code == 200:
        html_data = response.text

        # Got the page, so parse it, and try each selector
        fullsoup = BeautifulSoup(html_data, features="lxml")
        for selector in selectors:
            selected = fullsoup.select(selector)
            if selected:

                # Selected the article title, so assign it
                if len(selected) > 1:
                    logging.warning(
                        f"Selected more than one element using {selector} on soup from {citation_url}"
                    )
                title = selected[0].text
                try_wget = False
                break

    if try_wget:

        # Attempt to get the publication page using wget
        print(f"Trying wget")
        sleep(HTTPS_SLEEP)
        completed_process = subprocess.run(
            ["curl", "-L", citation_url], capture_output=True
        )
        html_data = completed_process.stdout

        # Got the page, so parse it, and search for the title
        fullsoup = BeautifulSoup(html_data, features="lxml")
        found = fullsoup.find_all("script")
        if found and len(found) > 4:
            m2 = p2.search(found[4].text)
            if m2:
                title = m2.group(1)

    print(f"Found title: '{title}' for citation URL: {citation_url}")

    return title


def get_pmid_for_title(title):
    """Search PubMed using a title to find the corresponding PMID.

    Parameters
    ----------
    title : str
       The title to use in the search

    Returns
    -------
    pmid : str
       The PubMed identifier found
    """
    # Need a default return value
    pmid = None

    # Search PubMed
    if title is None:
        return
    print(f"Getting PMID for title: '{title}'")
    search_url = EUTILS_URL + "esearch.fcgi"
    params = {
        "db": PUBMED,
        "term": title,
        "field": "title",
        "retmode": "json",
        # "retmax": 0,
        "email": NCBI_EMAIL,
        "api_key": NCBI_API_KEY,
    }
    sleep(NCBI_API_SLEEP)
    response = requests.get(search_url, params=parse.urlencode(params, safe=","))
    if response.status_code == 200:
        data = response.json()
        resultcount = int(data["esearchresult"]["count"])

        if resultcount > 1:
            # Response contains more than once result, so fetch each
            # PMID until title matches
            logging.warning(f"PubMed returned more than one result for title: {title}")
            for _pmid in data["esearchresult"]["idlist"]:
                _title = get_title_for_pmid(_pmid)
                if (
                    _title == title + "."
                ):  # PMID fetch includes period in title, title search does not
                    pmid = _pmid
                    break

        else:
            pmid = data["esearchresult"]["idlist"][0]

        print(f"Found PMID: {pmid} for title: '{title}'")

    elif response.status_code == 429:
        logging.error("Too many requests to NCBI API. Try again later, or use API key.")

    else:
        logging.error("Encountered error in searching PubMed: {response.status_code}")

    return pmid


def get_title_for_pmid(pmid):
    """Fetch from PubMed using a PMID to find the corresponding title.

    Parameters
    ----------
    pmid : str
       The PubMed identifier to use in the fetch

    Returns
    -------
    title : str
       The title fetched
    """
    # Need a default return value
    title = None

    # Fetch from PubMed
    fetch_url = EUTILS_URL + "efetch.fcgi"
    params = {
        "db": PUBMED,
        "id": pmid,
        "rettype": "xml",
        "email": NCBI_EMAIL,
        "api_key": NCBI_API_KEY,
    }
    sleep(NCBI_API_SLEEP)
    response = requests.get(fetch_url, params=parse.urlencode(params, safe=","))
    if response.status_code == 200:
        xml_data = response.text

        # Got the page, so parse it, and search for the title
        fullsoup = BeautifulSoup(xml_data, "xml")
        found = fullsoup.find("ArticleTitle")
        if found:
            title = found.text

    else:
        logging.error(
            f"Encountered error in fetching from PubMed: {response.status_code}"
        )

    return title


def get_and_download_dataset_h5ad_file(dataset_series):
    """Get the dataset filename and download the dataset file.

    Parameters
    ----------
    dataset_series : pd.Series
        A row from the dataset DataFrame

    Returns
    -------
    dataset : str
       The dataset filename
    """
    # Need a default return value
    dataset_filename = None

    # Get the dataset object
    collection_id = dataset_series.collection_id
    dataset_id = dataset_series.dataset_id
    dataset_url = f"{CELLXGENE_API_URL_BASE}/curation/v1/collections/{collection_id}/datasets/{dataset_id}"
    sleep(HTTPS_SLEEP)
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

    # Find H5AD files, if possible
    assets = data["assets"]
    for asset in assets:
        if asset["filetype"] != "H5AD":
            continue

        # Found an H5AD file, so download it, if needed
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


def run_nsforest(lung_datasets):
    """Run NSForest for each dataset file.

    Parameters
    ----------
    lung_datasets : pd.DataFrame
        DataFrame containing dataset descriptions with dataset
        filenames appended

    Returns
    -------
    None
    """
    # Run NSForest for each dataset file
    for dataset_h5ad_file in lung_datasets["dataset_h5ad_file"]:
        try:
            run_nsforest_on_file(dataset_h5ad_file)
        except Exception as ex:
            print(
                f"Could not run NSForest for unprocessed AnnData file: {dataset_h5ad_file}: {ex}"
            )


def run_nsforest_on_file(h5ad_filename, cluster_header="cell_type"):
    """Run NSForest using the specified dataset filename, and
    cluster_header.

    Parameters
    ----------
    h5ad_filename : str
       The dataset filename
    cluster_header : str
       The cluster header

    Returns
    -------
    None

    Notes
    -----
    - Some datasets have multiple annotations per sample
    (ex. "broad_cell_type" and "granular_cell_type"). NSForest can be
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
    # Assign results filename and directory
    pp_h5ad_filename = f"pp_{h5ad_filename}"
    results_dirname = h5ad_filename.split(".")[0]
    results_dirpath = f"{NSFOREST_DIR}/{results_dirname}"

    # Run NSForest if results do not exist
    if not os.path.exists(results_dirpath):
        os.makedirs(results_dirpath)

        print(f"Loading unprocessed AnnData file: {h5ad_filename}")
        h5ad_filepath = f"{CELLXGENE_DIR}/{h5ad_filename}"
        up_adata = sc.read_h5ad(h5ad_filepath)

        # TODO: Check validity of downsampling
        print("Calculating QC metrics")
        up_metrics = sc.pp.calculate_qc_metrics(up_adata)
        if up_metrics[1]["total_counts"].sum() > TOTAL_COUNTS:
            print("Downsampling unprocessed AnnData file")
            ds_adata = sc.pp.downsample_counts(
                up_adata, total_counts=TOTAL_COUNTS, copy=True
            )
        else:
            ds_adata = up_adata  # No need to copy

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

        pp_h5ad_filepath = f"{results_dirpath}/{pp_h5ad_filename}"
        print(f"Saving preprocessed AnnData file: {pp_h5ad_filepath}")
        pp_adata.write_h5ad(pp_h5ad_filepath)

        print(f"Running NSForest for preprocessed AnnData file: {pp_h5ad_filename}")
        results = nsforesting.NSForest(
            pp_adata,
            cluster_header,
            output_folder=f"{results_dirpath}/",
            outputfilename_prefix=cluster_header,
        )

    else:
        print(f"Completed NSForest for preprocessed AnnData file: {pp_h5ad_filename}")


def run_ontogpt(lung_datasets):
    """Run the OntoGPT pubmed-annotate function for each PMID
    associated with a dataset.

    Parameters
    ----------
    lung_datasets : pd.DataFrame
        DataFrame containing dataset descriptions with PMIDs appended

    Returns
    -------
    None
    """
    print("Running OntoGPT")
    with Pool(8) as p:
        p.map(run_ontogpt_pubmed_annotate, lung_datasets["citation_pmid"])


def run_ontogpt_pubmed_annotate(pmid):
    """Run the OntoGPT pubmed-annotate function for the specified PMID
    associated with a dataset.

    Parameters
    ----------
    pmid : str
       The PubMed identifier found

    Returns
    -------
    None
    """
    # Run OntoGPT pubmed-annotate function, if needed
    if pmid is None:
        return
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


def main():

    up_lung_obs, up_lung_datasets = get_lung_obs_and_datasets()
    pp_lung_datasets = append_titles_pmids_and_dataset_h5ad_files(up_lung_datasets)

    run_ontogpt(pp_lung_datasets)
    run_nsforest(pp_lung_datasets)

    return up_lung_obs, pp_lung_datasets


if __name__ == "__main__":
    __spec__ = None  # Workaround for Pool() in IPython
    up_lung_obs, pp_lung_datasets = main()
