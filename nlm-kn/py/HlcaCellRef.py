import argparse
import ast
from pathlib import Path
import random
import string

import ArangoDB as adb
import pandas as pd

from CellOntology import parse_term


ALPHABET = string.ascii_lowercase + string.digits


def get_uuid():
    return "".join(random.choices(ALPHABET, k=8))


def load_data(data_dir, data_fnm):
    data = pd.read_excel(data_dir / data_fnm, header=1, skiprows=0)
    data["uuid"] = [get_uuid() for idx in data.index]
    return data


def get_graph(db_name, graph_name):
    db = adb.create_or_get_database(db_name)
    adb_graph = adb.create_or_get_graph(db, graph_name)
    return adb_graph


def init_collections(adb_graph):

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


def add_vertices(vertex_collections):

    # Anatomic structures vertex
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


def add_vertices_and_edges_from_row(row, vertex_collections, edge_collections):

    anatomic_structure_vertex, publication_hlca_vertex, publication_cellref_vertex = (
        add_vertices(vertex_collections)
    )

    # Biomarker combination vertices
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

    # Cell set vertices
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

    # Cell type vertices
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
    # cell_type_cl_vertex = {}
    if not (
        pd.isna(row["CL_cell_type"])
        or pd.isna(row["CL_PURL"])
        # or pd.isna(row["Current CL definition"])
        # or pd.isna(row["Proposed addition to CL definition or annotation property."])
    ):
        if "http" in row["CL_PURL"]:
            _, number, _, _, _ = parse_term(row["CL_PURL"])

        elif ":" in row["CL_PURL"]:
            number = row["CL_PURL"].split(":")[1]
        _key = number
        cell_type_cl_vertex = vertex_collections["CL"].get(_key)
        if cell_type_cl_vertex:
            cell_type_cl_vertex["proposed definition"] = str(
                row["Proposed addition to CL definition or annotation property."]
            )
            vertex_collections["CL"].update(cell_type_cl_vertex)
        else:
            cell_type_cl_vertex = {}
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

    # Gene vertices
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

    # Biomarker combination IS_MARKER_FOR cell set triples
    triples = []
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

    # Cell type PART_OF anatomic structure triples
    triples.extend(
        [
            # (cell_type_hlca_vertex, {"name": "PART_OF"}, anatomic_structure_vertex),
            # (cell_type_cellref_vertex, {"name": "PART_OF"}, anatomic_structure_vertex),
            (cell_type_cl_vertex, {"name": "PART_OF"}, anatomic_structure_vertex),
        ]
    )

    # Cell set IS_INSTANCE cell type triples
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

    # Cell set EXPRESSES gene triples
    for gene in hlca_genes:
        _key = gene
        gene_vertex = vertex_collections["gene"].get(_key)
        triples.append((cell_set_hlca_vertex, {"name": "EXPRESSES"}, gene_vertex))
    for gene in cellref_genes:
        _key = gene
        gene_vertex = vertex_collections["gene"].get(_key)
        triples.append((cell_set_cellref_vertex, {"name": "EXPRESSES"}, gene_vertex))

    # Cell set SOURCE publication triples
    triples.extend(
        [
            (cell_set_hlca_vertex, {"name": "SOURCE"}, publication_hlca_vertex),
            (cell_set_cellref_vertex, {"name": "SOURCE"}, publication_cellref_vertex),
        ]
    )

    # Gene PART_OF biomarker combination triples
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

    for triple in triples:

        from_vertex = triple[0]  # subject
        predicate = triple[1]
        to_vertex = triple[2]  # object

        if from_vertex == {} or to_vertex == {}:
            continue

        edge = {
            "_key": from_vertex["_key"] + "-" + to_vertex["_key"],
            "_from": from_vertex["_id"],
            "_to": to_vertex["_id"],
        }
        edge.update(predicate)

        from_vertex_collection = from_vertex["_id"].split("/")[0]
        to_vertex_collection = to_vertex["_id"].split("/")[0]
        edge_name = from_vertex_collection + "-" + to_vertex_collection

        if not edge_collections[edge_name].has(edge):
            edge_collections[edge_name].insert(edge)


def main():

    parser = argparse.ArgumentParser(description="Load Manually Curated Data")
    parser.add_argument(
        "--data-dirname",
        default=Path("../data/nlm-kn"),
        help="Name of directory containing manually curated data",
    )
    parser.add_argument(
        "--data-filename",
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

    data = load_data(args.data_dirname, args.data_filename)

    adb_graph = get_graph(args.db_name, args.graph_name)

    vertex_collections, edge_collections = init_collections(adb_graph)

    for _, row in data.iterrows():

        add_vertices_and_edges_from_row(row, vertex_collections, edge_collections)
