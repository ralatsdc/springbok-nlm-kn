import logging
import os
from time import sleep
from traceback import print_exc
from urllib import parse

from bs4 import BeautifulSoup
import requests

DATA_DIR = "../data"

EUTILS_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
NCBI_EMAIL = os.environ.get("NCBI_EMAIL")
NCBI_API_KEY = os.environ.get("NCBI_API_KEY")
NCBI_API_SLEEP = 1
PUBMED = "pubmed"


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
        return pmid
    print(f"Getting PMID for title: '{title}'")
    search_url = EUTILS_URL + "esearch.fcgi"
    print(search_url)
    params = {
        "db": PUBMED,
        "term": title,
        "field": "title",
        "retmode": "json",
        # "retmax": 0,
        "email": NCBI_EMAIL,
        "api_key": NCBI_API_KEY,
    }
    print(params)
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
        logging.error(f"Encountered error in searching PubMed: {response.status_code}")

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
