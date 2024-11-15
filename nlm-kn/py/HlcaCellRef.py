import argparse
import ast
import json
from pathlib import Path
import random
import string

import gget
import ArangoDB as adb
import pandas as pd
import scanpy as sc

from CellOntology import parse_term


ALPHABET = string.ascii_lowercase + string.digits


def get_uuid():
    """Get an eight character random string.

    Parameters
    ----------
    None

    Returns
    -------
    An eight character random string.
    """
    return "".join(random.choices(ALPHABET, k=8))


def get_gene_name_to_ids_map():

    # TODO: Write documentation

    # Retrieve gene name to ids mapping
    return sc.queries.biomart_annotations(
        "hsapiens", ["ensembl_gene_id", "external_gene_name"], use_cache=True
    ).set_index("external_gene_name")


def map_gene_name_to_ids(name, nm2id):

    # TODO: Write documentation

    if name in nm2id.index:
        id = nm2id.loc[name, "ensembl_gene_id"]
        if isinstance(id, pd.core.series.Series):
            return id.to_list()
        else:
            return [id]
    else:
        print(f"Could not find ids for name: {name}")
        return []


def load_mdata(mdata_dirname, mdata_filename):
    """Load manually curated data.

    Parameters
    ----------
    mdata_dirname : str
        Name of directory containing manually curated data
    mdata_filename : str
        Name of Excel file containing manually curated data

    Returns
    -------
    mdata : pd.DataFrame
        DataFrame containing manually curated data
    """
    mdata = pd.read_excel(mdata_dirname / mdata_filename, header=1, skiprows=0)
    mdata["uuid"] = [get_uuid() for idx in mdata.index]
    return mdata


def init_gdata(mdata):

    # TODO: Write documentation

    gdata = {}
    gdata["cell_types"] = {}
    gdata["marker_names"] = set()
    gdata["marker_ids"] = set()

    nm2ids = get_gene_name_to_ids_map()

    for _, row in mdata.iterrows():

        if pd.isna(row["HLCA_cellset"]):
            continue
        cell_type = row["HLCA_cellset"]
        gdata["cell_types"][cell_type] = {}

        if pd.isna(row["HLCA_NSForestMarkers"]):
            continue
        marker_names = ast.literal_eval(row["HLCA_NSForestMarkers"])
        gdata["cell_types"][cell_type]["marker_names"] = {}
        gdata["marker_names"] |= set(marker_names)

        for marker_name in marker_names:
            marker_ids = map_gene_name_to_ids(marker_name, nm2ids)
            gdata["cell_types"][cell_type]["marker_names"][marker_name] = {}
            gdata["cell_types"][cell_type]["marker_names"][marker_name][
                "marker_ids"
            ] = marker_ids
            gdata["marker_ids"] |= set(marker_ids)

    gdata["marker_names"] = list(gdata["marker_names"])
    gdata["marker_ids"] = list(gdata["marker_ids"])

    return gdata


def load_gdata(gdata):

    # TODO: Write documentation

    resources = [
        "diseases",
        "drugs",
        "interactions",
        "tractability",
        "expression",
        "depmap",
    ]

    gdata["gene_ids"] = {}
    for id in gdata["marker_ids"][0:5]:
        gdata["gene_ids"][id] = {}

        gdata["gene_ids"][id]["resources"] = {}
        for resource in resources:
            try:
                gdata["gene_ids"][id]["resources"][resource] = gget.opentargets(
                    id, resource=resource, json=True, verbose=True
                )
            except Exception as exc:
                print(f"Could not gget resource: {resource} for id: {id}")
                gdata["gene_ids"][id]["resources"][resource] = {}


def get_graph(db_name, graph_name):
    """Get an graph from an ArangoDB database.

    Parameters
    ----------
    db_name : str
        Name of ArangoDB database
    graph_name : str
        Name of ArangoDB database graph

    Returns
    -------
    adb_graph : arango.graph.Graph
        ArangoDB database graph
    """
    db = adb.create_or_get_database(db_name)
    adb_graph = adb.create_or_get_graph(db, graph_name)
    return adb_graph


def init_collections(adb_graph):
    """Define and create vertex and edge collections.

    Parameters
    ----------
    adb_graph : arango.graph.Graph
        ArangoDB database graph

    Returns
    -------
    vertex_collections : dict
        A dictionary with vertex name keys containing
        arango.collection.VertexCollection instance values
    edge_collections : dict
        A dictionary with edge name keys containing
        arango.collection.EdgeCollection instance values
    """
    # Define and create vertex collections
    vertex_collections = {}
    vertex_names = [
        "anatomic_structure",
        "biomarker_combination",
        "cell_set",
        "CL",
        "gene",
        "publication",
    ]
    for vertex_name in vertex_names:
        collection = adb.create_or_get_vertex_collection(adb_graph, vertex_name)
        vertex_collections[vertex_name] = collection

    # Define and create edge collections
    edge_collections = {}
    vertex_name_pairs = [
        ("biomarker_combination", "cell_set"),
        ("CL", "anatomic_structure"),
        ("cell_set", "CL"),
        ("cell_set", "gene"),
        ("cell_set", "publication"),
        ("CL", "CL"),
        ("gene", "biomarker_combination"),
    ]
    for vertex_name_pair in vertex_name_pairs:
        from_vertex = vertex_name_pair[0]
        to_vertex = vertex_name_pair[1]
        collection, edge_name = adb.create_or_get_edge_collection(
            adb_graph, from_vertex, to_vertex
        )
        edge_collections[edge_name] = collection

    return vertex_collections, edge_collections


def insert_vertices(vertex_collections):
    """Add the anatomic structure and publication vertices.

    Parameters
    ----------
    vertex_collections : dict
        A dictionary with vertex name keys containing
        arango.collection.VertexCollection instance values

    Returns
    -------
    anatomic_structure_vertex : dict
        The ArangoDB anatomic structure vertex document
    publication_hlca_vertex : dict
        The ArangoDB HLCA publication vertex document
    publication_cellref_vertex,
        The ArangoDB CellRef publication vertex document
    """
    # Anatomic structure vertex
    _key = "Organ_12345"
    if not vertex_collections["anatomic_structure"].has(_key):
        vertex_collections["anatomic_structure"].insert({"_key": _key, "name": "lung"})
    anatomic_structure_vertex = vertex_collections["anatomic_structure"].get(_key)

    # Publication vertices
    _key = "HLCA_2023_Sikkema"
    if not vertex_collections["publication"].has(_key):
        vertex_collections["publication"].insert(
            {
                "_key": _key,
                "name": "HLCA_2023_Sikkema_doi.org/10.1038/s41591-023-02327-2",
            }
        )
    publication_hlca_vertex = vertex_collections["publication"].get(_key)
    _key = "cellRef_2023_Guo"
    if not vertex_collections["publication"].has(_key):
        vertex_collections["publication"].insert(
            {
                "_key": _key,
                "name": "cellRef_2023_Guo_doi.org/10.1038/s41467-023-40173-5",
            }
        )
    publication_cellref_vertex = vertex_collections["publication"].get(_key)

    return (
        anatomic_structure_vertex,
        publication_hlca_vertex,
        publication_cellref_vertex,
    )


def insert_vertices_and_edges_from_row(row, vertex_collections, edge_collections):
    """Add vertices and edges from the given row of the DataFrame containing manually curated data.

    Parameters
    ----------
    row : pd.Series
        Row Series from the DataFrame containing manually curated data

    Returns
    -------
    None
    """
    # Add the anatomic structure and publication vertices
    anatomic_structure_vertex, publication_hlca_vertex, publication_cellref_vertex = (
        insert_vertices(vertex_collections)
    )

    # Add biomarker combination vertices
    hlca_genes = []
    biomarker_combination_hlca_vertex = {}
    cellref_genes = []
    biomarker_combination_cellref_vertex = {}
    if row["HLCA_NSForestMarkers"] == row["CellRef_NSForestMarkers"] and not pd.isna(
        row["HLCA_NSForestMarkers"]
    ):
        hlca_genes = ast.literal_eval(row["HLCA_NSForestMarkers"])
        cellref_genes = hlca_genes
        _key = "hlca-cellref-" + row["uuid"]
        if not vertex_collections["biomarker_combination"].get(_key):
            vertex_collections["biomarker_combination"].insert(
                {
                    "_key": _key,
                    "name": hlca_genes,
                }
            )
        biomarker_combination_hlca_vertex = vertex_collections[
            "biomarker_combination"
        ].get(_key)
        biomarker_combination_cellref_vertex = biomarker_combination_hlca_vertex

    else:
        if not pd.isna(row["HLCA_NSForestMarkers"]):
            hlca_genes = ast.literal_eval(row["HLCA_NSForestMarkers"])
            _key = "hlca-" + row["uuid"]
            if not vertex_collections["biomarker_combination"].has(_key):
                vertex_collections["biomarker_combination"].insert(
                    {"_key": _key, "name": hlca_genes}
                )
            biomarker_combination_hlca_vertex = vertex_collections[
                "biomarker_combination"
            ].get(_key)

        if not pd.isna(row["CellRef_NSForestMarkers"]):
            cellref_genes = ast.literal_eval(row["CellRef_NSForestMarkers"])
            _key = "cellref-" + row["uuid"]
            if not vertex_collections["biomarker_combination"].has(_key):
                vertex_collections["biomarker_combination"].insert(
                    {
                        "_key": _key,
                        "name": cellref_genes,
                    }
                )
            biomarker_combination_cellref_vertex = vertex_collections[
                "biomarker_combination"
            ].get(_key)

    # Add cell set vertices
    cell_set_hlca_vertex = {}
    if not pd.isna(row["HLCA_cellset"]):
        _key = "hlca-" + row["uuid"]
        if not vertex_collections["cell_set"].has(_key):
            vertex_collections["cell_set"].insert(
                {
                    "_key": _key,
                    "name": row["HLCA_cellset"],
                }
            )
        cell_set_hlca_vertex = vertex_collections["cell_set"].get(_key)
    cell_set_cellref_vertex = {}
    if not pd.isna(row["cellref_cellset"]):
        _key = "cellref-" + row["uuid"]
        if not vertex_collections["cell_set"].has(_key):
            vertex_collections["cell_set"].insert(
                {
                    "_key": _key,
                    "name": row["cellref_cellset"],
                }
            )
        cell_set_cellref_vertex = vertex_collections["cell_set"].get(_key)

    # Add cell type vertices
    # TODO: Decide whether we need these cell types, or not
    # cell_type_hlca_vertex = {}
    # if not pd.isna(row["Cell_type_HLCA"]):
    #     _key = "hlca-" + row["uuid"]
    #     if not vertex_collections["CL"].has(_key):
    #         vertex_collections["CL"].insert(
    #             {
    #                 "_key": _key,
    #                 "name": row["Cell_type_HLCA"],
    #             }
    #         )
    #     cell_type_hlca_vertex = vertex_collections["CL"].get(_key)
    # cell_type_cellref_vertex = {}
    # if not pd.isna(row["Cell_type_cellref"]):
    #     _key = "cellref-" + row["uuid"]
    #     if not vertex_collections["CL"].has(_key):
    #         vertex_collections["CL"].insert(
    #             {
    #                 "_key": _key,
    #                 "name": row["Cell_type_cellref"],
    #             }
    #         )
    #     cell_type_cellref_vertex = vertex_collections["CL"].get(_key)
    cell_type_cl_vertex = {}
    if not (
        pd.isna(row["CL_cell_type"])
        or pd.isna(row["CL_PURL"])
        # TODO: Decide if we need to test if these definitions are present
        # or pd.isna(row["Current CL definition"])
        # or pd.isna(row["Proposed addition to CL definition or annotation property."])
    ):
        # Get the number from the current CL term
        if "http" in row["CL_PURL"]:
            _, number, _, _, _ = parse_term(row["CL_PURL"])

        elif ":" in row["CL_PURL"]:
            number = row["CL_PURL"].split(":")[1]

        # Update the CL term vertex, if it exists
        _key = number
        cell_type_cl_vertex = vertex_collections["CL"].get(_key)
        if cell_type_cl_vertex:
            cell_type_cl_vertex["proposed definition"] = str(
                row["Proposed addition to CL definition or annotation property."]
            )
            vertex_collections["CL"].update(cell_type_cl_vertex)
        else:
            cell_type_cl_vertex = {}
        # TODO: Remove
        # if not vertex_collections["CL"].has(_key):
        #     vertex_collections["CL"].insert(
        #         {
        #             "_key": _key,
        #             "name": row["CL_cell_type"],
        #             "purl": row["CL_PURL"],
        #             "current definition": str(row["Current CL definition"]),
        #             "proposed definition": str(
        #                 row["Proposed addition to CL definition or annotation property."]
        #             ),
        #         }
        #     )
        # cell_type_cl_vertex = vertex_collections["CL"].get(_key)

    # Add gene vertices
    for gene in hlca_genes + cellref_genes:
        _key = gene
        if not vertex_collections["gene"].has(_key):
            vertex_collections["gene"].insert(
                {
                    "_key": _key,
                    "name": gene,
                }
            )
            gene_vertex = vertex_collections["gene"].get(_key)

    # Initialize triples
    triples = []

    # Define and collect (biomarker combination, IS_MARKER_FOR, cell set) triples
    triples.extend(
        [
            (
                biomarker_combination_hlca_vertex,
                {"name": "IS_MARKER_FOR", "Fbeta": row["HLCA_Fbeta"]},
                cell_set_hlca_vertex,
            ),
            (
                biomarker_combination_cellref_vertex,
                {"name": "IS_MARKER_FOR", "Fbeta": row["CellRef_Fbeta"]},
                cell_set_cellref_vertex,
            ),
        ]
    )

    # Define and collect (cell type, PART_OF, anatomic structure) triples
    triples.extend(
        [
            # (cell_type_hlca_vertex, {"name": "PART_OF"}, anatomic_structure_vertex),
            # (cell_type_cellref_vertex, {"name": "PART_OF"}, anatomic_structure_vertex),
            (cell_type_cl_vertex, {"name": "PART_OF"}, anatomic_structure_vertex),
        ]
    )

    # Define and collect (cell set, IS_INSTANCE, cell type) triples
    if (
        row["predicate_HLCA"] == "skos:exactMatch"
        or row["predicate_HLCA"] == "skos:relatedMatch"
    ):
        triples.extend(
            [
                # (cell_set_hlca_vertex, {"name": "IS_INSTANCE"}, cell_type_hlca_vertex),
                (
                    cell_set_hlca_vertex,
                    {"name": "IS_INSTANCE", "match": row["predicate_HLCA"]},
                    cell_type_cl_vertex,
                ),
                # (cell_set_cellref_vertex, {"name": "IS_INSTANCE"}, cell_type_cellref_vertex),
                (
                    cell_set_cellref_vertex,
                    {"name": "IS_INSTANCE", "match": row["predicate_HLCA"]},
                    cell_type_cl_vertex,
                ),
            ]
        )

    # Define and collect (cell set, EXPRESSES, gene) triples
    for gene in hlca_genes:
        _key = gene
        gene_vertex = vertex_collections["gene"].get(_key)
        triples.append((cell_set_hlca_vertex, {"name": "EXPRESSES"}, gene_vertex))
    for gene in cellref_genes:
        _key = gene
        gene_vertex = vertex_collections["gene"].get(_key)
        triples.append((cell_set_cellref_vertex, {"name": "EXPRESSES"}, gene_vertex))

    # Define and collect (cell set, SOURCE, publication) triples
    triples.extend(
        [
            (cell_set_hlca_vertex, {"name": "SOURCE"}, publication_hlca_vertex),
            (cell_set_cellref_vertex, {"name": "SOURCE"}, publication_cellref_vertex),
        ]
    )

    # Define and collect (gene, PART_OF, biomarker combination) triples
    for gene in hlca_genes:
        _key = gene
        gene_vertex = vertex_collections["gene"].get(_key)
        triples.append(
            (gene_vertex, {"name": "PART_OF"}, biomarker_combination_hlca_vertex)
        )
    for gene in cellref_genes:
        _key = gene
        gene_vertex = vertex_collections["gene"].get(_key)
        triples.append(
            (
                gene_vertex,
                {"name": "PART_OF"},
                biomarker_combination_cellref_vertex,
            )
        )

    # Insert all triples
    for triple in triples:

        from_vertex = triple[0]  # subject
        predicate = triple[1]
        to_vertex = triple[2]  # object

        # All vertices are valid or an empty dictionary
        if from_vertex == {} or to_vertex == {}:
            continue

        # Define edge, noting that all predicates are arbitrary dictionaries
        edge = {
            "_key": from_vertex["_key"] + "-" + to_vertex["_key"],
            "_from": from_vertex["_id"],
            "_to": to_vertex["_id"],
        }
        edge.update(predicate)

        # Insert triple
        from_vertex_collection = from_vertex["_id"].split("/")[0]
        to_vertex_collection = to_vertex["_id"].split("/")[0]
        edge_name = from_vertex_collection + "-" + to_vertex_collection
        if not edge_collections[edge_name].has(edge):
            edge_collections[edge_name].insert(edge)


def main():

    parser = argparse.ArgumentParser(description="Load Manually Curated Data")
    parser.add_argument(
        "--mdata-dirname",
        default=Path("../data/nlm-kn"),
        help="Name of directory containing manually curated data",
    )
    parser.add_argument(
        "--mdata-filename",
        default="HLCA_CellRef_matching_ver3_import1.xlsm",
        help="Name of file containing manually curated data",
    )
    group = parser.add_argument_group(
        "Cell Ontology (CL)", "Assume this version of the CL has been loaded"
    )
    exclusive_group = group.add_mutually_exclusive_group(required=True)
    exclusive_group.add_argument(
        "--slim", action="store_true", help="Assume the slim ontology has been loaded"
    )
    exclusive_group.add_argument(
        "--full", action="store_true", help="Assume the full ontology has been loaded"
    )
    group.add_argument(
        "--include-bnodes",
        action="store_true",
        help="Assume BNodes were included when loading",
    )

    args = parser.parse_args()

    if args.slim:
        db_name = "BioPortal-Slim"
        graph_name = "CL-Slim"

    if args.full:
        db_name = "BioPortal-Full"
        graph_name = "CL-Full"

    if args.include_bnodes:
        db_name += "-BNodes"
        graph_name += "-BNodes"

    print(f"Loading {args.mdata_dirname / args.mdata_filename}")
    mdata = load_mdata(args.mdata_dirname, args.mdata_filename)

    gdata = init_gdata(mdata)

    load_gdata(gdata)

    gdata_filename = args.mdata_filename.replace(".xlsm", ".json")
    with open(args.mdata_dirname / gdata_filename, "w") as ofp:
        json.dump(gdata, ofp, indent=4)

    return mdata, gdata

    # print(f"Getting graph {graph_name} from {db_name}")
    # adb_graph = get_graph(db_name, graph_name)

    # print("Defining and creating vertex and edge collections")
    # vertex_collections, edge_collections = init_collections(adb_graph)

    # print("Inserting vertices and edges form each row of the manually curated data")
    # for _, row in mdata.iterrows():
    #     insert_vertices_and_edges_from_row(row, vertex_collections, edge_collections)


if __name__ == "__main__":
    mdata, gdata = main()
