import os
import subprocess
from traceback import print_exc

DATA_DIR = "../data"

ONTOGPT_DIR = f"{DATA_DIR}/ontogpt"


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
