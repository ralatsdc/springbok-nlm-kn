import ast
from glob import glob
import json
import logging
import os
import re
import yaml

from arango import ArangoClient
import pandas as pd


ARANGO_URL = "http://localhost:8529"
ARANGO_CLIENT = ArangoClient(hosts=ARANGO_URL)
SYS_DB = ARANGO_CLIENT.db("_system", username="root", password="")

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

    # Create NCBI-Cell database, if needed
    if not SYS_DB.has_database(NCBI_CELL_DB):
        SYS_DB.create_database(NCBI_CELL_DB)

    # Connect to NCBI-Cell database
    db = ARANGO_CLIENT.db(NCBI_CELL_DB, username="root", password="")

    return db


def delete_database():

    # Delete NCBI-Cell database, if needed
    if SYS_DB.has_database(NCBI_CELL_DB):
        SYS_DB.delete_database(NCBI_CELL_DB)


def create_or_get_graph(db):

    # Create, or get the graph
    if db.has_graph(NCBI_CELL_GRAPH):
        graph = db.graph(NCBI_CELL_GRAPH)
    else:
        graph = db.create_graph(NCBI_CELL_GRAPH)

    return graph, db


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

    return cellxgene, graph


def create_and_populate_or_get_vertex_collection_nsforest(graph):

    # Create, or get the nsforest vertex collection
    if graph.has_vertex_collection(NSFOREST_COLLECTION):
        nsforest = graph.vertex_collection(NSFOREST_COLLECTION)
    else:
        nsforest = graph.create_vertex_collection(NSFOREST_COLLECTION)

    # Read each NSForest results file
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

    return nsforest, graph


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

    return ontogpt, graph


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
            cell_id = m.group(1)
            gpt["cell_id"] = cell_id

            d = {k: gpt[k] for k in ("id", "citation_pmid")}
            d["_key"] = cell_id
            if not cell.has(d["_key"]):
                cell.insert(d)

        else:
            gpt["cell_id"] = None
        ontogpt.update(gpt)

    return cell, graph, nsforest, ontogpt


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

    return gene, graph, nsforest


def create_or_get_edge_collection_cellxgene_cell(graph):

    if not graph.has_edge_definition(f"{CELLXGENE_COLLECTION}-{CELL_COLLECTION}"):
        cellxgene_cell = graph.create_edge_definition(
            edge_collection=f"{CELLXGENE_COLLECTION}-{CELL_COLLECTION}",
            from_vertex_collections=[f"{CELLXGENE_COLLECTION}"],
            to_vertex_collections=[f"{CELL_COLLECTION}"],
        )
    else:
        cellxgene_cell = graph.edge_collection(
            f"{CELLXGENE_COLLECTION}-{CELL_COLLECTION}"
        )

    return cellxgene_cell, graph


def create_or_get_edge_collection_ontogpt_cell(graph):

    if not graph.has_edge_definition(f"{ONTOGPT_COLLECTION}-{CELL_COLLECTION}"):
        ontogpt_cell = graph.create_edge_definition(
            edge_collection=f"{ONTOGPT_COLLECTION}-{CELL_COLLECTION}",
            from_vertex_collections=[f"{ONTOGPT_COLLECTION}"],
            to_vertex_collections=[f"{CELL_COLLECTION}"],
        )
    else:
        ontogpt_cell = graph.edge_collection(f"{ONTOGPT_COLLECTION}-{CELL_COLLECTION}")

    return ontogpt_cell, graph


def create_or_get_edge_collection_cell_gene(graph):

    if not graph.has_edge_definition(f"{CELL_COLLECTION}-{GENE_COLLECTION}"):
        cell_gene = graph.create_edge_definition(
            edge_collection=f"{CELL_COLLECTION}-{GENE_COLLECTION}",
            from_vertex_collections=[f"{CELL_COLLECTION}"],
            to_vertex_collections=[f"{GENE_COLLECTION}"],
        )
    else:
        cell_gene = graph.edge_collection(f"{CELL_COLLECTION}-{GENE_COLLECTION}")

    return cell_gene, graph


def insert_cellxgene_cell_edges(cellxgene, cell, cellxgene_cell):

    for cxg in cellxgene.all():
        cxg_key = cxg["_key"]
        print(
            f"Finding edges to {CELL_COLLECTION} from {CELLXGENE_COLLECTION} document with key: {cxg_key}"
        )

        found = False
        for cll in cell.all():
            if "dataset_id" in cll and cll["dataset_id"] == cxg["dataset_id"]:
                found = True
                cll_key = cll["_key"]
                doc = {
                    "_key": f"{cxg_key}-{cll_key}",
                    "_from": f"{CELLXGENE_COLLECTION}/{cxg_key}",
                    "_to": f"{CELL_COLLECTION}/{cll_key}",
                }
                if not cellxgene_cell.has(doc):
                    print(
                        f"Inserting edge to {CELL_COLLECTION} document with key: {cll_key} from {CELLXGENE_COLLECTION} document with key: {cxg_key}"
                    )
                    cellxgene_cell.insert(doc)

        if not found:
            print(
                f"No edges to {CELL_COLLECTION} from {CELLXGENE_COLLECTION} document with key: {cxg_key}"
            )

    return cellxgene, cell, cellxgene_cell


def insert_ontogpt_cell_edges(ontogpt, cell, ontogpt_cell):

    for gpt in ontogpt.all():
        gpt_key = gpt["_key"]
        print(
            f"Finding edges to {CELL_COLLECTION} from {ONTOGPT_COLLECTION} document with key: {gpt_key}"
        )

        found = False
        cll_key = gpt["cell_id"]
        if cll_key is not None and cell.has(cll_key):
            found = True

            doc = {
                "_key": f"{gpt_key}-{cll_key}",
                "_from": f"{ONTOGPT_COLLECTION}/{gpt_key}",
                "_to": f"{CELL_COLLECTION}/{cll_key}",
            }
            if not ontogpt_cell.has(doc):
                print(
                    f"Inserting edge to {CELL_COLLECTION} document with key: {cll_key} from {ONTOGPT_COLLECTION} document with key: {gpt_key}"
                )
                ontogpt_cell.insert(doc)

        if not found:
            print(
                f"No edges to {CELL_COLLECTION} from {ONTOGPT_COLLECTION} document with key: {gpt_key}"
            )

    return ontogpt, cell, ontogpt_cell


def insert_cell_gene_edges(nsforest, cell_gene):

    for nsf in nsforest.all():
        cll_key = nsf["clusterName"]

        for gn_key in ast.literal_eval(nsf["NSForest_markers"]):
            doc = {
                "_key": f"{cll_key}-{gn_key}",
                "_from": f"{CELL_COLLECTION}/{cll_key}",
                "_to": f"{GENE_COLLECTION}/{gn_key}",
            }
            if not cell_gene.has(doc):
                print(
                    f"Inserting edge to {GENE_COLLECTION} document with key: {gn_key} from {CELL_COLLECTION} document with key: {cll_key}"
                )
                cell_gene.insert(doc)

    return nsforest, cell_gene


def main():
    pass


if __name__ == "__main__":

    delete_database()
    db = create_or_get_database()

    graph, db = create_or_get_graph(db)

    cellxgene, graph = create_and_populate_or_get_vertex_collection_cellxgene(graph)
    nsforest, graph = create_and_populate_or_get_vertex_collection_nsforest(graph)
    ontogpt, graph = create_and_populate_or_get_vertex_collection_ontogpt(graph)
    cell, graph, nsforest, ontogpt = create_and_populate_or_get_vertex_collection_cell(
        graph, nsforest, ontogpt
    )
    gene, graph, nsforest = create_and_populate_or_get_vertex_collection_gene(
        graph, nsforest
    )

    cellxgene_cell, graph = create_or_get_edge_collection_cellxgene_cell(graph)
    ontogpt_cell, graph = create_or_get_edge_collection_ontogpt_cell(graph)
    cell_gene, graph = create_or_get_edge_collection_cell_gene(graph)

    cellxgene, cell, cellxgene_cell = insert_cellxgene_cell_edges(
        cellxgene, cell, cellxgene_cell
    )
    ontogpt, cell, ontogpt_cell = insert_ontogpt_cell_edges(ontogpt, cell, ontogpt_cell)
    nsforest, cell_gene = insert_cell_gene_edges(nsforest, cell_gene)
