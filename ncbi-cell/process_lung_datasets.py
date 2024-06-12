# https://chanzuckerberg.github.io/cellxgene-census/notebooks/analysis_demo/comp_bio_explore_and_load_lung_data.html
import logging
from multiprocessing.pool import Pool, ThreadPool
import os
import re
import requests
import subprocess
from time import sleep
from urllib import parse

from bs4 import BeautifulSoup
import cellxgene_census
import pandas as pd


EUTILS_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
EMAIL = "raymond.leclair@gmail.com"
NCBI_API_KEY = os.environ.get("NCBI_API_KEY")
NCBI_API_SLEEP = 1
PUBMED = "pubmed"
PUBMEDCENTRAL = "pmc"


def get_lung_obs_and_datasets():

    lung_obs_parquet = "lung_obs.parquet"
    lung_datasets_parquet = "lung_datasets.parquet"

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

    citations = [c for c in lung_datasets.citation]
    with ThreadPool(8) as p:
        titles = p.map(get_title, citations)

    # for citation in citations:
    #     title = get_title(citation)
    #     titles.append(title)

    return list(set([title for title in titles if title is not None]))


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

    with Pool(8) as p:
        pmids = p.map(get_pmid_for_title, titles)

    # for title in titles:
    #     pmid = get_pmid_for_title(title)
    #     pmids.append(pmid)

    return list(set([pmid for pmid in pmids if pmid is not None]))


def run_ontogpt_pubmed_annotate(pmid):
    """
    run_ontogpt_pubmed_annotate("38540357")
    """
    if not os.path.exists(f"{pmid}.out"):
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
                f"{pmid}.out",
            ],
        )

        print(f"Completed ontogpt pubmed-annotate for PMID: {pmid}")

    else:
        print(f"Ontogpt pubmed-annotate out for PMID: {pmid} exists")


def run_ontogpt(pmids):

    with ThreadPool(8) as p:
        p.map(run_ontogpt_pubmed_annotate, pmids)


if __name__ == "__main__":

    lung_obs, lung_datasets = get_lung_obs_and_datasets()

    titles = get_titles(lung_datasets)

    pmids = get_pmids(titles)

    run_ontogpt(pmids)
