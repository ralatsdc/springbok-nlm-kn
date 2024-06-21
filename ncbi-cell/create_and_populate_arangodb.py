import ast
from glob import glob
import json
import os
import re
import yaml

from arango import ArangoClient
import pandas as pd


ARANGO_URL = "http://localhost:8529"
ARANGO_CLIENT = ArangoClient(hosts=ARANGO_URL)

DATA_DIR = "data"

NCBI_CELL_DIR = f"{DATA_DIR}/ncbi-cell"
NCBI_CELL_DB = "ncbi-cell"
NCBI_CELL_GRAPH = "ncbi-cell"

CELLXGENE_COLLECTION = "cellxgene"

NSFOREST_DIR = f"{DATA_DIR}/nsforest"
NSFOREST_COLLECTION = "nsforest"

ONTOGPT_DIR = f"{DATA_DIR}/ontogpt"
ONTOGPT_COLLECTION = "ontogpt"

CELL_COLLECTION = "cell"
GENE_COLLECTION = "gene"


def create_or_get_database():

    # Connect to "_system" database as root user.
    sys_db = ARANGO_CLIENT.db("_system", username="root", password="")

    # Create NCBI-Cell database, if needed
    if not sys_db.has_database(NCBI_CELL_DB):
        sys_db.create_database(NCBI_CELL_DB)

    # Connect to NCBI-Cell database
    db = ARANGO_CLIENT.db(NCBI_CELL_DB, username="root", password="")

    return db


def delete_database():

    # Connect to "_system" database as root user.
    sys_db = ARANGO_CLIENT.db("_system", username="root", password="")

    # Delete NCBI-Cell database, if needed
    if sys_db.has_database(NCBI_CELL_DB):
        sys_db.delete_database(NCBI_CELL_DB)


def create_or_get_graph(db):

    # Create, or get the graph
    if db.has_graph(NCBI_CELL_GRAPH):
        graph = db.graph(NCBI_CELL_GRAPH)
    else:
        graph = db.create_graph(NCBI_CELL_GRAPH)

    return graph


def create_and_populate_or_get_vertex_collection_cellxgene(graph):

    # Create or get the cellxgne vertex collection
    if graph.has_vertex_collection(CELLXGENE_COLLECTION):
        cellxgene = graph.vertex_collection(CELLXGENE_COLLECTION)
    else:
        cellxgene = graph.create_vertex_collection(CELLXGENE_COLLECTION)

    # Read preprocessed lung datasets parquet
    pp_lung_datasets_parquet = f"{NCBI_CELL_DIR}/pp_lung_datasets.parquet"
    pp_lung_datasets = pd.read_parquet(pp_lung_datasets_parquet)

    # Insert documents using collection_id as _key
    pp_lung_datasets["_key"] = pp_lung_datasets["collection_id"]
    for index, row in pp_lung_datasets.iterrows():
        if not cellxgene.has(row["_key"]):
            cellxgene.insert(json.loads(row.to_json()))

    return cellxgene


def create_and_populate_or_get_vertex_collection_nsforest(graph):

    # Create, or get the nsforest vertex collection
    if graph.has_vertex_collection(NSFOREST_COLLECTION):
        nsforest = graph.vertex_collection(NSFOREST_COLLECTION)
    else:
        nsforest = graph.create_vertex_collection(NSFOREST_COLLECTION)

    # Read each OntoGPT results file
    for fn in glob(f"{NSFOREST_DIR}/*/*.csv"):
        df = pd.read_csv(fn)

        # Append the dataset_id
        dataset_id = os.path.basename(os.path.dirname(fn))
        df["dataset_id"] = dataset_id

        # Insert documents using clusterName as _key
        df["_key"] = df["clusterName"]
        for index, row in df.iterrows():
            if not nsforest.has(row["_key"]):
                nsforest.insert(json.loads(row.to_json()))

    return nsforest


def create_and_populate_or_get_vertex_collection_ontogpt(graph):

    # Create, or get the ontogpt vertex collection
    if graph.has_vertex_collection(ONTOGPT_COLLECTION):
        ontogpt = graph.vertex_collection(ONTOGPT_COLLECTION)
    else:
        ontogpt = graph.create_vertex_collection(ONTOGPT_COLLECTION)

    # Load each OntoGPT output file
    for fn in glob(f"{ONTOGPT_DIR}/*.out"):
        with open(fn, "r") as f:
            y = yaml.safe_load(f)

        # Append the PMID
        pmid, _ = os.path.splitext(os.path.basename(fn))
        y["extracted_object"]["citation_pmid"] = pmid

        # Use the extracted_object dictionary as the document, then
        # insert it using its id as _key
        d = y["extracted_object"]
        d["_key"] = d["id"]
        if not ontogpt.has(d["_key"]):
            ontogpt.insert(d)

    return ontogpt


def create_and_populate_or_get_vertex_collection_cell(graph, nsforest, ontogpt):

    # Create, or get the cell vertex collection
    if graph.has_vertex_collection(CELL_COLLECTION):
        cell = graph.vertex_collection(CELL_COLLECTION)
    else:
        cell = graph.create_vertex_collection(CELL_COLLECTION)

    # Subset each nsforest document, then insert it using its
    # clusterName as _key
    for nsf in nsforest.all():
        d = {k: nsf[k] for k in ("clusterName", "dataset_id")}
        d["_key"] = d["clusterName"]
        if not cell.has(d["_key"]):
            cell.insert(d)

    # Search for cells in each ontogpt document, and if found, subset
    # the document, then insert it using the first match
    p = re.compile("'(CL:\d*)'")
    for gpt in ontogpt.all():
        m = p.search(str(gpt))
        if m:
            d = {k: gpt[k] for k in ("id", "citation_pmid")}
            d["_key"] = m.group(1)
            if not cell.has(d["_key"]):
                cell.insert(d)

    return cell


def create_and_populate_or_get_vertex_collection_gene(graph, nsforest):

    # Create, or get the gene vertex collection
    if graph.has_vertex_collection(GENE_COLLECTION):
        gene = graph.vertex_collection(GENE_COLLECTION)
    else:
        gene = graph.create_vertex_collection(GENE_COLLECTION)

    # Subset each document, then insert it using each of
    # NSForest_markers as _key
    for nsf in nsforest.all():
        d = {k: nsf[k] for k in ("clusterName", "dataset_id")}
        for mrk in ast.literal_eval(nsf["NSForest_markers"]):
            d["_key"] = mrk
            if not gene.has(d["_key"]):
                gene.insert(d)

    return gene


def abc(cellxgene, nsforest, ontogpt):

    for cxg in cellxgene.all():

        for nsf in nsforest.all():
            if nsf["dataset_id"] == cxg["dataset_id"]:
                break

        for gpt in ontogpt.all():
            if gpt["citation_pmid"] == cxg["citation_pmid"]:
                break

        break

    return


def main():

    delete_database()
    db = create_or_get_database()

    graph = create_or_get_graph(db)

    cellxgene = create_and_populate_or_get_vertex_collection_cellxgene(graph)
    nsforest = create_and_populate_or_get_vertex_collection_nsforest(graph)
    ontogpt = create_and_populate_or_get_vertex_collection_ontogpt(graph)
    cell = create_and_populate_or_get_vertex_collection_cell(graph, nsforest, ontogpt)
    gene = create_and_populate_or_get_vertex_collection_gene(graph, nsforest)

    return cellxgene, nsforest, ontogpt, cell, gene


if __name__ == "__main__":
    cellxgene, nsforest, ontogpt, cell, gene = main()
