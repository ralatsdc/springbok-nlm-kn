import ast
from glob import glob
from io import StringIO
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
NCBI_CELL_DB = "ncbi-cell-2024-06-27"
NCBI_CELL_GRAPH = "ncbi-cell"

CELLXGENE_COLLECTION = "cellxgene"

NSFOREST_DIR = f"{DATA_DIR}/nsforest-2024-06-27"
NSFOREST_COLLECTION = "nsforest"

ONTOGPT_DIR = f"{DATA_DIR}/ontogpt"
ONTOGPT_COLLECTION = "ontogpt"

CELL_COLLECTION = "cell"
GENE_COLLECTION = "gene"


def create_or_get_database():
    """Create or get NCBI-Cell database.

    Parameters
    ----------
    None

    Returns
    -------
    db : arango.database.StandardDatabase
        NCBI-Cell database
    """
    # Create NCBI-Cell database, if needed
    if not SYS_DB.has_database(NCBI_CELL_DB):
        SYS_DB.create_database(NCBI_CELL_DB)

    # Connect to NCBI-Cell database
    db = ARANGO_CLIENT.db(NCBI_CELL_DB, username="root", password="")

    return db


def delete_database():
    """Delete NCBI-Cell database.

    Parameters
    ----------
    None

    Returns
    -------
    None
    """
    # Delete NCBI-Cell database, if needed
    if SYS_DB.has_database(NCBI_CELL_DB):
        SYS_DB.delete_database(NCBI_CELL_DB)


def create_or_get_graph(db):
    """Create or get NCBI-Cell database graph.

    Parameters
    ----------
    db : arango.database.StandardDatabase
        NCBI-Cell database

    Returns
    -------
    graph : arango.graph.Graph
        NCBI-Cell database graph
    db : arango.database.StandardDatabase
        NCBI-Cell database
    """
    # Create, or get the graph
    if db.has_graph(NCBI_CELL_GRAPH):
        graph = db.graph(NCBI_CELL_GRAPH)
    else:
        graph = db.create_graph(NCBI_CELL_GRAPH)

    return graph, db


def delete_graph(db):
    """Delete NCBI-Cell database graph.

    Parameters
    ----------
    db : arango.database.StandardDatabase
        NCBI-Cell database

    Returns
    -------
    """
    # Delete the graph
    if db.has_graph(NCBI_CELL_GRAPH):
        db.delete_graph(NCBI_CELL_GRAPH)


def create_and_populate_or_get_vertex_collection_cellxgene(graph):
    """Create or get the CELLXGENE vertex collection, read the
    preprocessed datasets parquet, and insert a vertex for each row of
    the datasets DataFrame.

    Parameters
    ----------
    graph : arango.graph.Graph
        NCBI-Cell database graph

    Returns
    -------
    cellxgene : arango.collection.VertexCollection
        CELLXGENE vertex collection
    graph : arango.graph.Graph
        NCBI-Cell database graph
    """
    # Create or get the cellxgene vertex collection
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
    """Create or get the NSForest vertex collection, read each
    NSForest results file, and insert a vertex for results DataFrame.

    Parameters
    ----------
    graph : arango.graph.Graph
        NCBI-Cell database graph

    Returns
    -------
    nsforest : arango.collection.VertexCollection
        NSForest vertex collection
    graph : arango.graph.Graph
        NCBI-Cell database graph
    """
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

        # Insert documents using dataset_id as _key
        _key = dataset_id
        if not nsforest.has(_key):
            d = dataframe_to_doc(df, _key)
            nsforest.insert(d)

    return nsforest, graph


def create_and_populate_or_get_vertex_collection_ontogpt(graph):
    """Create or get the OntoGPT vertex collection, read each OntoGPT
    output file, and insert a vertex for output object.

    Parameters
    ----------
    graph : arango.graph.Graph
        NCBI-Cell database graph

    Returns
    -------
    ontogpt : arango.collection.VertexCollection
        OntoGPT vertex collection
    graph : arango.graph.Graph
        NCBI-Cell database graph
    """
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


def dataframe_to_doc(df, _key):
    """Convert a Pandas DataFrame to an ArangoDB document using the
    specified _key.

    Parameters
    ----------
    df : pd.DataFrame
        Pandas DataFrame

    Returns
    -------
    doc : dict
        ArangoDB document with key needed by ArangoDB
    """
    doc = json.loads(df.to_json())
    doc["_key"] = _key
    return doc


def doc_to_dataframe(doc):
    """Convert ArangoDB document to Pandas DataFrame.

    Parameters
    ----------
    doc : dict
        ArangoDB document

    Returns
    -------
    df : pd.DataFrame
        ArangoDB document without keys added by ArangoDB as Pandas
        DataFrame
    """
    for key in ["_key", "_id", "_rev"]:
        if key in doc:
            doc.pop(key)
    df = pd.read_json(StringIO(json.dumps(doc)))
    return df


def create_and_populate_or_get_vertex_collection_cell(graph, nsforest, ontogpt):
    """Create or get the cell vertex collection, and insert a vertex
    corresponding to each cell identified in the nsforest and ontogpt
    vertex collections.

    Parameters
    ----------
    graph : arango.graph.Graph
        NCBI-Cell database graph
    nsforest : arango.collection.VertexCollection
        NSForest vertex collection
    ontogpt : arango.collection.VertexCollection
        OntoGPT vertex collection

    Returns
    -------
    cell : arango.collection.VertexCollection
        cell vertex collection
    graph : arango.graph.Graph
        NCBI-Cell database graph
    nsforest : arango.collection.VertexCollection
        NSForest vertex collection
    ontogpt : arango.collection.VertexCollection
        OntoGPT vertex collection
    """
    # Create, or get the cell vertex collection
    if graph.has_vertex_collection(CELL_COLLECTION):
        cell = graph.vertex_collection(CELL_COLLECTION)
    else:
        cell = graph.create_vertex_collection(CELL_COLLECTION)

    # Read each nsforest vertex JSON to create a DataFrame
    for nsf in nsforest.all():
        df = doc_to_dataframe(nsf)

        # Consider each row of the DataFrame
        for index, row in df.iterrows():
            _key = row["clusterName"].replace(" ", "-").replace(",", ":")

            # Insert or update a cell vertex using the row clusterName
            # as key, collecting all dataset_ids corresponding to the
            # cell vertex
            if not cell.has(_key):
                d = {
                    "_key": _key,
                    "clusterName": row["clusterName"],
                    "dataset_ids": [row["dataset_id"]],
                }
                cell.insert(d)

            else:
                d = cell.get(_key)
                d["dataset_ids"].append(row["dataset_id"])
                cell.update(d)

    # Search for cell ontology identifiers in each ontogpt document,
    # and if found, subset the document, then insert it using the
    # first match
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

        # Update the ontgpt vertex with the cell_id
        ontogpt.update(gpt)

    return cell, graph, nsforest, ontogpt


def create_and_populate_or_get_vertex_collection_gene(graph, nsforest):
    """Create or get the gene vertex collection, and insert a vertex
    corresponding to each marker for each cell identified in the
    nsforest vertex collection.

    Parameters
    ----------
    graph : arango.graph.Graph
        NCBI-Cell database graph
    nsforest : arango.collection.VertexCollection
        NSForest vertex collection

    Returns
    -------
    graph : arango.graph.Graph
        NCBI-Cell database graph
    nsforest : arango.collection.VertexCollection
        NSForest vertex collection
    """
    # Create, or get the gene vertex collection
    if graph.has_vertex_collection(GENE_COLLECTION):
        gene = graph.vertex_collection(GENE_COLLECTION)
    else:
        gene = graph.create_vertex_collection(GENE_COLLECTION)

    # Read each nsforest vertex JSON to create a DataFrame
    for nsf in nsforest.all():
        df = doc_to_dataframe(nsf)

        # Consider each row of the DataFrame
        for index, row in df.iterrows():

            # Consider each marker in the row
            for mrk in ast.literal_eval(row["NSForest_markers"]):
                _key = mrk

                # Insert or update a gene vertex using the marker as
                # key, collecting all clusterNames and dataset_ids
                # corresponding to the gene vertex
                if not gene.has(_key):
                    d = {
                        "_key": _key,
                        "clusterNames": [row["clusterName"]],
                        "dataset_ids": [row["dataset_id"]],
                    }
                    gene.insert(d)

                else:
                    d = gene.get(_key)
                    d["clusterNames"].append(row["clusterName"])
                    d["dataset_ids"].append(row["dataset_id"])
                    gene.update(d)

    return gene, graph, nsforest


def create_or_get_edge_collection_cellxgene_cell(graph):
    """Create, or get the from CELLXGENE to cell edge collection.

    Parameters
    ----------
    graph : arango.graph.Graph
        NCBI-Cell database graph

    Returns
    -------
    cellxgene_cell : arango.collection.EdgeCollection
        From CELLXGENE to cell edge collection
    graph : arango.graph.Graph
        NCBI-Cell database graph
    """
    # Create, or get the from cellxgene to cell edge collection
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
    """Create, or get the from OntoGPT to cell edge collection.

    Parameters
    ----------
    graph : arango.graph.Graph
        NCBI-Cell database graph

    Returns
    -------
    ontogpt_cell : arango.collection.EdgeCollection
        From OntoGPT to cell edge collection
    graph : arango.graph.Graph
        NCBI-Cell database graph
    """
    # Create, or get the from ontogpt to cell edge collection
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
    """Create, or get the from cell to gene edge collection.

    Parameters
    ----------
    graph : arango.graph.Graph
        NCBI-Cell database graph

    Returns
    -------
    cell_gene : arango.collection.EdgeCollection
        From cell to gene edge collection
    graph : arango.graph.Graph
        NCBI-Cell database graph
    """
    # Create, or get the from cell to gene edge collection
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
    """Insert an edge from any CELLXGENE vertex to a cell vertex if
    they share their dataset_id.

    Parameters
    ----------
    cellxgene : arango.collection.VertexCollection
        CELLXGENE vertex collection
    cell : arango.collection.VertexCollection
        cell vertex collection
    cellxgene_cell : arango.collection.EdgeCollection
        From CELLXGENE to cell edge collection

    Returns
    -------
    cellxgene : arango.collection.VertexCollection
        CELLXGENE vertex collection
    cell : arango.collection.VertexCollection
        cell vertex collection
    cellxgene_cell : arango.collection.EdgeCollection
        From CELLXGENE to cell edge collection
    """
    # Consider each cellxgene vertex
    for cxg in cellxgene.all():
        cxg_key = cxg["_key"]  # collection_id
        print(
            f"Finding edges to {CELL_COLLECTION} from {CELLXGENE_COLLECTION} document with key: {cxg_key}"
        )

        # Consider each cell vertex
        found = False
        for cll in cell.all():

            # Insert an edge from the cellxgene vertex to the cell
            # vertex if they share their dataset_id
            if "dataset_ids" in cll and cxg["dataset_id"] in cll["dataset_ids"]:
                found = True
                cll_key = cll["_key"]  # clusterName with replacements
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
    """Insert an edge from an OntoGPT vertex to a cell vertex if the
    OntoGPT vertex cell_id is the cell key. Currently, this cannot
    happen.

    Parameters
    ----------
    ontogpt : arango.collection.VertexCollection
        OntoGPT vertex collection
    cell : arango.collection.VertexCollection
        cell vertex collection
    ontogpt_cell : arango.collection.EdgeCollection
        From OntoGPT to cell edge collection

    Returns
    -------
    ontogpt : arango.collection.VertexCollection
        OntoGPT vertex collection
    cell : arango.collection.VertexCollection
        cell vertex collection
    ontogpt_cell : arango.collection.EdgeCollection
        From OntoGPT to cell edge collection
    """
    # Consider each ontogpt vertex
    for gpt in ontogpt.all():
        gpt_key = gpt["_key"]  # id
        print(
            f"Finding edges to {CELL_COLLECTION} from {ONTOGPT_COLLECTION} document with key: {gpt_key}"
        )

        # Insert an edge from an ontogpt vertex to a cell vertex if
        # the ontogpt vertex cell_id is the cell key. Currently, this
        # cannot happen.
        cll_key = gpt["cell_id"]  # cell ontology identifier
        if cll_key is not None and cell.has(cll_key):
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

        else:
            print(
                f"No edges to {CELL_COLLECTION} from {ONTOGPT_COLLECTION} document with key: {gpt_key}"
            )

    return ontogpt, cell, ontogpt_cell


def insert_cell_gene_edges(nsforest, cell_gene):
    """Use the NSForest vertex collection to insert an edge from every
    cell vertex to the marker gene vertex.

    Parameters
    ----------
    nsforest : arango.collection.VertexCollection
        NSForest vertex collection
    cell_gene : arango.collection.EdgeCollection
        From cell to gene edge collection

    Returns
    ----------
    nsforest : arango.collection.VertexCollection
        NSForest vertex collection
    cell_gene : arango.collection.EdgeCollection
        From cell to gene edge collection
    """
    # Read each nsforest vertex JSON to create a DataFrame
    for nsf in nsforest.all():
        df = doc_to_dataframe(nsf)

        # Consider each row of the DataFrame which corresponds to a
        # cell vertex
        for index, row in df.iterrows():
            cll_key = row["clusterName"].replace(" ", "-").replace(",", ":")

            # Consider each marker in the row which corresponds to a
            # gene vertex
            for gn_key in ast.literal_eval(row["NSForest_markers"]):

                # Insert an edge from the cell vertex to the gene
                # vertex, if needed
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
    # TODO: Use this function to clean up variables in scope
    pass


if __name__ == "__main__":

    # Create a new database
    delete_database()
    db = create_or_get_database()

    # Create a graph
    graph, db = create_or_get_graph(db)

    # Create and populate vertex collections
    cellxgene, graph = create_and_populate_or_get_vertex_collection_cellxgene(graph)
    nsforest, graph = create_and_populate_or_get_vertex_collection_nsforest(graph)
    ontogpt, graph = create_and_populate_or_get_vertex_collection_ontogpt(graph)
    cell, graph, nsforest, ontogpt = create_and_populate_or_get_vertex_collection_cell(
        graph, nsforest, ontogpt
    )
    gene, graph, nsforest = create_and_populate_or_get_vertex_collection_gene(
        graph, nsforest
    )

    # Create edge collections
    cellxgene_cell, graph = create_or_get_edge_collection_cellxgene_cell(graph)
    ontogpt_cell, graph = create_or_get_edge_collection_ontogpt_cell(graph)
    cell_gene, graph = create_or_get_edge_collection_cell_gene(graph)

    # Insert edges
    cellxgene, cell, cellxgene_cell = insert_cellxgene_cell_edges(
        cellxgene, cell, cellxgene_cell
    )
    ontogpt, cell, ontogpt_cell = insert_ontogpt_cell_edges(ontogpt, cell, ontogpt_cell)
    nsforest, cell_gene = insert_cell_gene_edges(nsforest, cell_gene)
