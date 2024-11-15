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
    """Query BioMart to map gene names to ids.

    Parameters
    ----------
    None

    Returns
    -------
    pd.DataFrame
        DataFrame indexed by gene name containing gene id
    """
    return sc.queries.biomart_annotations(
        "hsapiens", ["ensembl_gene_id", "external_gene_name"], use_cache=True
    ).set_index("external_gene_name")


def map_gene_name_to_ids(name, nm2id):
    """Map a gene name to a gene id list.

    Parameters
    ----------
    name : str
        Gene name
    nm2id : pd.DataFrame
        DataFrame indexed by name containing id

    Returns
    -------
    list
        Gene ids
    """
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
    mdata_dirname : Path
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
    """Initialize gget data.

    Parameters
    ----------
    mdata : pd.DataFrame
        DataFrame containing manually curated data

    Returns
    -------
    gdata : dict
        Dictionary containing HLCA cell types, marker names, marker
        ids for each marker name, and all unique marker ids
    """
    nm2ids = get_gene_name_to_ids_map()
    gdata = {}
    gdata["cell_types"] = {}
    gdata["marker_ids"] = set()
    for _, row in mdata.iterrows():

        # Collect HLCA cell types
        if pd.isna(row["HLCA_cellset"]):
            continue
        cell_type = row["HLCA_cellset"]
        gdata["cell_types"][cell_type] = {}

        # Collect marker names
        if pd.isna(row["HLCA_NSForestMarkers"]):
            continue
        marker_names = ast.literal_eval(row["HLCA_NSForestMarkers"])
        gdata["cell_types"][cell_type]["marker_names"] = {}
        for marker_name in marker_names:

            # Collect marker id
            marker_ids = map_gene_name_to_ids(marker_name, nm2ids)
            gdata["cell_types"][cell_type]["marker_names"][marker_name] = {}
            gdata["cell_types"][cell_type]["marker_names"][marker_name][
                "marker_ids"
            ] = marker_ids
            gdata["marker_ids"] |= set(marker_ids)

    # Marker ids are more convenient as a list
    gdata["marker_ids"] = list(gdata["marker_ids"])

    return gdata


def load_gdata(mdata, gdata_dirname, gdata_filename):
    """Call gget for each gene id, collecting all specified resources.

    Parameters
    ----------
    mdata : pd.DataFrame
        DataFrame containing manually curated data
    gdata_dirname : Path
        Name of directory containing gget data
    gdata_filename : str
        Name of JSON file containing gget data

    Returns
    -------
    gdata : dict
        Dictionary containing HLCA cell types, marker names, marker
        ids for each marker name, all specified resources for each
        marker id, and all unique marker ids
    """
    # Create and dump, or load gget data
    if not (gdata_dirname / gdata_filename).exists():
        # Create gget data
        gdata = init_gdata(mdata)
        resources = [
            "diseases",
            "drugs",
            "interactions",
            "tractability",
            "expression",
            "depmap",
        ]

        # Collect gene ids
        gdata["gene_ids"] = {}
        for id in gdata["marker_ids"]:
            gdata["gene_ids"][id] = {}

            # Collect resources
            gdata["gene_ids"][id]["resources"] = {}
            for resource in resources:
                try:
                    gdata["gene_ids"][id]["resources"][resource] = gget.opentargets(
                        id, resource=resource, json=True, verbose=True
                    )
                except Exception as exc:
                    print(f"Could not gget resource: {resource} for id: {id}")
                    gdata["gene_ids"][id]["resources"][resource] = {}

        # Dump gget data
        with open(gdata_dirname / gdata_filename, "w") as fp:
            json.dump(gdata, fp, indent=4)

    else:

        # Load gget data
        with open(gdata_dirname / gdata_filename, "r") as fp:
            gdata = json.load(fp)

    return gdata


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
        "gene_name",
        "publication",
        "gene_id",
        "disease",
        "drug_product",
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
        ("cell_set", "gene_name"),
        ("cell_set", "publication"),
        ("CL", "CL"),
        ("gene_name", "biomarker_combination"),
        ("gene_name", "gene_id"),
        ("gene_id", "disease"),
        ("drug_product", "gene_id"),
        ("drug_product", "disease"),
    ]
    for vertex_name_pair in vertex_name_pairs:
        from_vertex = vertex_name_pair[0]
        to_vertex = vertex_name_pair[1]
        collection, edge_name = adb.create_or_get_edge_collection(
            adb_graph, from_vertex, to_vertex
        )
        edge_collections[edge_name] = collection

    return vertex_collections, edge_collections


def insert_anatomic_structure_and_publication_vertices(vertex_collections):
    """Insert the single anatomic structure and two publication
    vertices.

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
        vertex_collections["anatomic_structure"].insert({"_key": _key, "label": "lung"})
    anatomic_structure_vertex = vertex_collections["anatomic_structure"].get(_key)

    # Publication vertices
    _key = "HLCA_2023_Sikkema"
    if not vertex_collections["publication"].has(_key):
        vertex_collections["publication"].insert(
            {
                "_key": _key,
                "label": "HLCA_2023_Sikkema_doi.org/10.1038/s41591-023-02327-2",
            }
        )
    publication_hlca_vertex = vertex_collections["publication"].get(_key)
    _key = "cellRef_2023_Guo"
    if not vertex_collections["publication"].has(_key):
        vertex_collections["publication"].insert(
            {
                "_key": _key,
                "label": "cellRef_2023_Guo_doi.org/10.1038/s41467-023-40173-5",
            }
        )
    publication_cellref_vertex = vertex_collections["publication"].get(_key)

    return (
        anatomic_structure_vertex,
        publication_hlca_vertex,
        publication_cellref_vertex,
    )


def insert_triples(triples, edge_collections):
    """Insert triples into the ArangoDB graph.

    Parameters
    ----------
    triples : list(tuple)
    edge_collections : dict
        A dictionary with edge name keys containing
        arango.collection.EdgeCollection instance values

    Returns
    -------
    None
    """
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


def insert_vertices_and_edges_from_mdata_row(row, vertex_collections, edge_collections):
    """Insert vertices and edges from each row of the manually curated
    data.

    Parameters
    ----------
    row : pd.Series
        Row Series from the DataFrame containing manually curated data

    Returns
    -------
    None
    """
    # Insert the anatomic structure and publication vertices
    anatomic_structure_vertex, publication_hlca_vertex, publication_cellref_vertex = (
        insert_anatomic_structure_and_publication_vertices(vertex_collections)
    )

    # Insert biomarker combination vertices
    hlca_gene_names = []
    biomarker_combination_hlca_vertex = {}
    cellref_gene_names = []
    biomarker_combination_cellref_vertex = {}
    if row["HLCA_NSForestMarkers"] == row["CellRef_NSForestMarkers"] and not pd.isna(
        row["HLCA_NSForestMarkers"]
    ):
        hlca_gene_names = ast.literal_eval(row["HLCA_NSForestMarkers"])
        cellref_gene_names = hlca_gene_names
        _key = "hlca-cellref-" + row["uuid"]
        if not vertex_collections["biomarker_combination"].get(_key):
            vertex_collections["biomarker_combination"].insert(
                {
                    "_key": _key,
                    "label": hlca_gene_names,
                }
            )
        biomarker_combination_hlca_vertex = vertex_collections[
            "biomarker_combination"
        ].get(_key)
        biomarker_combination_cellref_vertex = biomarker_combination_hlca_vertex

    else:
        if not pd.isna(row["HLCA_NSForestMarkers"]):
            hlca_gene_names = ast.literal_eval(row["HLCA_NSForestMarkers"])
            _key = "hlca-" + row["uuid"]
            if not vertex_collections["biomarker_combination"].has(_key):
                vertex_collections["biomarker_combination"].insert(
                    {"_key": _key, "label": hlca_gene_names}
                )
            biomarker_combination_hlca_vertex = vertex_collections[
                "biomarker_combination"
            ].get(_key)

        if not pd.isna(row["CellRef_NSForestMarkers"]):
            cellref_gene_names = ast.literal_eval(row["CellRef_NSForestMarkers"])
            _key = "cellref-" + row["uuid"]
            if not vertex_collections["biomarker_combination"].has(_key):
                vertex_collections["biomarker_combination"].insert(
                    {
                        "_key": _key,
                        "label": cellref_gene_names,
                    }
                )
            biomarker_combination_cellref_vertex = vertex_collections[
                "biomarker_combination"
            ].get(_key)

    # Insert cell set vertices
    cell_set_hlca_vertex = {}
    if not pd.isna(row["HLCA_cellset"]):
        _key = "hlca-" + row["uuid"]
        if not vertex_collections["cell_set"].has(_key):
            vertex_collections["cell_set"].insert(
                {
                    "_key": _key,
                    "label": row["HLCA_cellset"],
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
                    "label": row["cellref_cellset"],
                }
            )
        cell_set_cellref_vertex = vertex_collections["cell_set"].get(_key)

    # Insert cell type vertices
    # TODO: Decide whether we need these cell types, or not
    # cell_type_hlca_vertex = {}
    # if not pd.isna(row["Cell_type_HLCA"]):
    #     _key = "hlca-" + row["uuid"]
    #     if not vertex_collections["CL"].has(_key):
    #         vertex_collections["CL"].insert(
    #             {
    #                 "_key": _key,
    #                 "label": row["Cell_type_HLCA"],
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
    #                 "label": row["Cell_type_cellref"],
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
        #             "label": row["CL_cell_type"],
        #             "purl": row["CL_PURL"],
        #             "current definition": str(row["Current CL definition"]),
        #             "proposed definition": str(
        #                 row["Proposed addition to CL definition or annotation property."]
        #             ),
        #         }
        #     )
        # cell_type_cl_vertex = vertex_collections["CL"].get(_key)

    # Insert gene vertices
    for gene_name in hlca_gene_names + cellref_gene_names:
        _key = gene_name
        if not vertex_collections["gene_name"].has(_key):
            vertex_collections["gene_name"].insert(
                {
                    "_key": _key,
                    "label": gene_name,
                }
            )
            gene_name_vertex = vertex_collections["gene_name"].get(_key)

    # Initialize triples
    triples = []

    # Define and collect (biomarker combination, IS_MARKER_FOR, cell set) triples
    triples.extend(
        [
            (
                biomarker_combination_hlca_vertex,
                {"label": "IS_MARKER_FOR", "Fbeta": row["HLCA_Fbeta"]},
                cell_set_hlca_vertex,
            ),
            (
                biomarker_combination_cellref_vertex,
                {"label": "IS_MARKER_FOR", "Fbeta": row["CellRef_Fbeta"]},
                cell_set_cellref_vertex,
            ),
        ]
    )

    # Define and collect (cell type, PART_OF, anatomic structure) triples
    triples.extend(
        [
            # (cell_type_hlca_vertex, {"label": "PART_OF"}, anatomic_structure_vertex),
            # (cell_type_cellref_vertex, {"label": "PART_OF"}, anatomic_structure_vertex),
            (cell_type_cl_vertex, {"label": "PART_OF"}, anatomic_structure_vertex),
        ]
    )

    # Define and collect (cell set, IS_INSTANCE, cell type) triples
    if (
        row["predicate_HLCA"] == "skos:exactMatch"
        or row["predicate_HLCA"] == "skos:relatedMatch"
    ):
        triples.extend(
            [
                # (cell_set_hlca_vertex, {"label": "IS_INSTANCE"}, cell_type_hlca_vertex),
                (
                    cell_set_hlca_vertex,
                    {"label": "IS_INSTANCE", "match": row["predicate_HLCA"]},
                    cell_type_cl_vertex,
                ),
                # (cell_set_cellref_vertex, {"label": "IS_INSTANCE"}, cell_type_cellref_vertex),
                (
                    cell_set_cellref_vertex,
                    {"label": "IS_INSTANCE", "match": row["predicate_HLCA"]},
                    cell_type_cl_vertex,
                ),
            ]
        )

    # Define and collect (cell set, EXPRESSES, gene_name) triples
    for gene_name in hlca_gene_names:
        _key = gene_name
        gene_name_vertex = vertex_collections["gene_name"].get(_key)
        triples.append((cell_set_hlca_vertex, {"label": "EXPRESSES"}, gene_name_vertex))
    for gene_name in cellref_gene_names:
        _key = gene_name
        gene_name_vertex = vertex_collections["gene_name"].get(_key)
        triples.append(
            (cell_set_cellref_vertex, {"label": "EXPRESSES"}, gene_name_vertex)
        )

    # Define and collect (cell set, SOURCE, publication) triples
    triples.extend(
        [
            (cell_set_hlca_vertex, {"label": "SOURCE"}, publication_hlca_vertex),
            (cell_set_cellref_vertex, {"label": "SOURCE"}, publication_cellref_vertex),
        ]
    )

    # Define and collect (gene_name, PART_OF, biomarker combination) triples
    for gene_name in hlca_gene_names:
        _key = gene_name
        gene_name_vertex = vertex_collections["gene_name"].get(_key)
        triples.append(
            (gene_name_vertex, {"label": "PART_OF"}, biomarker_combination_hlca_vertex)
        )
    for gene_name in cellref_gene_names:
        _key = gene_name
        gene_name_vertex = vertex_collections["gene_name"].get(_key)
        triples.append(
            (
                gene_name_vertex,
                {"label": "PART_OF"},
                biomarker_combination_cellref_vertex,
            )
        )

    # Insert all triples
    insert_triples(triples, edge_collections)


def insert_vertices_and_edges_from_gdata_for_gene_name_and_id(
    gdata, gene_name, gene_id, vertex_collections, edge_collections
):
    """Insert vertices and edges for each gene name and id from the
    gget data.

    Parameters
    ----------
    gdata
    gene_name
    gene_id
    vertex_collections
    edge_collections

    Returns
    -------
    None
    """
    # Initialize triples
    triples = []

    # Get gene name vertex
    gene_name_vertex = {}
    if vertex_collections["gene_name"].has(gene_name):
        gene_name_vertex = vertex_collections["gene_name"].get(gene_name)

    # Insert gene id vertex
    _key = gene_id
    if not vertex_collections["gene_id"].has(_key):
        vertex_collections["gene_id"].insert(
            {
                "_key": _key,
                "label": gene_id,
            }
        )
    gene_id_vertex = vertex_collections["gene_id"].get(_key)

    # Define and collect (gene name, HAS, gene id) triples
    triples.append((gene_name_vertex, {"label": "HAS"}, gene_id_vertex))

    # Get diseases and drugs
    diseases = []
    drugs = []
    if gene_id in gdata["gene_ids"]:
        diseases = gdata["gene_ids"][gene_id]["resources"]["diseases"]
        drugs = gdata["gene_ids"][gene_id]["resources"]["drugs"]

    for disease in diseases:

        # Insert disease vertex
        _key = disease["id"]
        if not vertex_collections["disease"].has(_key):
            disease["_key"] = _key
            disease["label"] = disease["name"]
            vertex_collections["disease"].insert(disease)
        disease_vertex = vertex_collections["disease"].get(_key)

        # Define and collect (gene id, IS_BASIS_FOR_CONDITION, disease) triples
        triples.append(
            (
                gene_id_vertex,
                {"label": "IS_BASIS_FOR_CONDITION"},
                disease_vertex,
            )
        )

    for drug in drugs:

        # Insert drug product vertex
        _key = drug["id"]
        if not vertex_collections["drug_product"].has(_key):
            drug["_key"] = _key
            drug["label"] = drug["name"]
            vertex_collections["drug_product"].insert(drug)
        drug_product_vertex = vertex_collections["drug_product"].get(_key)

        # Define and collect (drug product, MOLECULARLY_INTERACTS_WITH, gene id) triples
        triples.append(
            (
                drug_product_vertex,
                {"label": "MOLECULARLY_INTERACTS_WITH"},
                gene_id_vertex,
            )
        )

        # Get diseases vertex
        _key = drug["disease_id"]
        if vertex_collections["disease"].has(_key):
            disease_vertex = vertex_collections["disease"].get(_key)

            # Define and collect (drug product, IS_SUBSTANCE_THAT_TREATS, disease) triples
            triples.append(
                (
                    drug_product_vertex,
                    {"label": "IS_SUBSTANCE_THAT_TREATS"},
                    disease_vertex,
                )
            )

    # Insert all triples
    insert_triples(triples, edge_collections)


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

    gdata_dirname = args.mdata_dirname
    gdata_filename = args.mdata_filename.replace(".xlsm", ".json")

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

    print(f"Loading {gdata_dirname / gdata_filename}")
    gdata = load_gdata(mdata, gdata_dirname, gdata_filename)

    print(f"Getting graph {graph_name} from {db_name}")
    adb_graph = get_graph(db_name, graph_name)

    print("Defining and creating vertex and edge collections")
    vertex_collections, edge_collections = init_collections(adb_graph)

    n_rows = mdata.shape[0]
    for i_row, row in mdata.iterrows():
        print("Inserting vertices and edges from row {i_orw} (of {n_row})of the manually curated data")
        insert_vertices_and_edges_from_mdata_row(
            row, vertex_collections, edge_collections
        )

    print("Inserting vertices and edges for each gene name and id from the gget data")
    for _, d in gdata["cell_types"].items():
        for gene_name, e in d["marker_names"].items():
            for gene_id in e["marker_ids"]:
                insert_vertices_and_edges_from_gdata_for_gene_name_and_id(
                    gdata, gene_name, gene_id, vertex_collections, edge_collections
                )


if __name__ == "__main__":
    main()
