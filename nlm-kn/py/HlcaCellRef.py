import ast
from pathlib import Path
import random
import string

import ArangoDB as adb
import pandas as pd


alphabet = string.ascii_lowercase + string.digits
def get_uuid():
    return "".join(random.choices(alphabet, k=8))


data_dir = Path("/Users/raymondleclair/Projects/NLM/NLM-KN/cell2Knowledge/data")
data_fnm = "HLCA_CellRef_matching_ver3_import1.xlsm"
data = pd.read_excel(data_dir / data_fnm, header=1, skiprows=0)
data["uuid"] = [get_uuid() for idx in data.index]


db_name = "hlca_cellref"
graph_name = "hlca_cellref"
adb.delete_database(db_name)
db = adb.create_or_get_database(db_name)
adb_graph = adb.create_or_get_graph(db, graph_name)


vertex_collections = {}
vertex_names = [
    "anatomic_structure",
    "biomarker_combination",
    "cell_set",
    "cell_type",
    "gene",
    "publication",
]
for vertex_name in vertex_names:
    collection = adb.create_or_get_vertex_collection(adb_graph, vertex_name)
    vertex_collections[vertex_name] = collection


edge_collections = {}
vertex_name_pairs = [
    ("biomarker_combination", "cell_set"),
    ("cell_type", "anatomic_structure"),
    ("cell_set", "cell_type"),
    ("cell_set", "gene"),
    ("cell_set", "publication"),
    ("cell_type", "cell_type"),
    ("gene", "biomarker_combination"),
]
for vertex_name_pair in vertex_name_pairs:
    from_vertex = vertex_name_pair[0]
    to_vertex = vertex_name_pair[1]
    collection, edge_name = adb.create_or_get_edge_collection(
        adb_graph, from_vertex, to_vertex
    )
    edge_collections[edge_name] = collection


# Anatomic structures vertex
_key = "Organ_12345"
vertex_collections["anatomic_structure"].insert({"_key": _key, "name": "lung"})
anatomic_structure_vertex = vertex_collections["anatomic_structure"].get(_key)


# Publications vertices
_key = "HLCA_2023_Sikkema"
vertex_collections["publication"].insert(
    {
        "_key": _key,
        "name": "HLCA_2023_Sikkema_doi.org/10.1038/s41591-023-02327-2",
    }
)
publication_hlca_vertex = vertex_collections["publication"].get(_key)
_key = "cellRef_2023_Guo"
vertex_collections["publication"].insert(
    {
        "_key": _key,
        "name": "cellRef_2023_Guo_doi.org/10.1038/s41467-023-40173-5",
    }
)
publication_cellref_vertex = vertex_collections["publication"].get(_key)


for _, row in data.iterrows():

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
            vertex_collections["biomarker_combination"].insert(
                {"_key": _key, "name": hlca_genes}
            )
            biomarker_combination_hlca_vertex = vertex_collections[
                "biomarker_combination"
            ].get(_key)

        if not pd.isna(row["CellRef_NSForestMarkers"]):
            cellref_genes = ast.literal_eval(row["CellRef_NSForestMarkers"])
            _key = "cellref-" + row["uuid"]
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
        vertex_collections["cell_set"].insert(
            {
                "_key": _key,
                "name": row["cellref_cellset"],
            }
        )
        cell_set_cellref_vertex = vertex_collections["cell_set"].get(_key)

    # Cell type vertices
    cell_type_hlca_vertex = {}
    if not pd.isna(row["Cell_type_HLCA"]):
        _key = "hlca-" + row["uuid"]
        vertex_collections["cell_type"].insert(
            {
                "_key": _key,
                "name": row["Cell_type_HLCA"],
            }
        )
        cell_type_hlca_vertex = vertex_collections["cell_type"].get(_key)
    cell_type_cellref_vertex = {}
    if not pd.isna(row["Cell_type_cellref"]):
        _key = "cellref-" + row["uuid"]
        vertex_collections["cell_type"].insert(
            {
                "_key": _key,
                "name": row["Cell_type_cellref"],
            }
        )
        cell_type_cellref_vertex = vertex_collections["cell_type"].get(_key)
    cell_type_cl_vertex = {}
    if not (
        pd.isna(row["CL_cell_type"])
        or pd.isna(row["CL_PURL"])
        or pd.isna(row["Current CL definition"])
        or pd.isna(row["Proposed addition to CL definition or annotation property."])
    ):
        _key = "cl-" + row["uuid"]
        vertex_collections["cell_type"].insert(
            {
                "_key": _key,
                "name": row["CL_cell_type"],
                "purl": row["CL_PURL"],
                "current definition": row["Current CL definition"],
                "proposed definition": row[
                    "Proposed addition to CL definition or annotation property."
                ],
            }
        )
        cell_type_cl_vertex = vertex_collections["cell_type"].get(_key)

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
            (biomarker_combination_hlca_vertex, "IS_MARKER_FOR", cell_set_hlca_vertex),
            (
                biomarker_combination_cellref_vertex,
                "IS_MARKER_FOR",
                cell_set_cellref_vertex,
            ),
        ]
    )

    # Cell type PART_OF anatomic structure triples
    triples.extend(
        [
            (cell_type_hlca_vertex, "PART_OF", anatomic_structure_vertex),
            (cell_type_cellref_vertex, "PART_OF", anatomic_structure_vertex),
        ]
    )

    # Cell set IS_INSTANCE cell type triples
    triples.extend(
        [
            (cell_set_hlca_vertex, "IS_INSTANCE", cell_type_hlca_vertex),
            (cell_set_hlca_vertex, "IS_INSTANCE", cell_type_cl_vertex),
            (cell_set_cellref_vertex, "IS_INSTANCE", cell_type_cellref_vertex),
            (cell_set_cellref_vertex, "IS_INSTANCE", cell_type_cl_vertex),
        ]
    )

    # Cell set EXPRESSES gene triples
    for gene in hlca_genes:
        gene_vertex = vertex_collections["gene"].get(gene)
        triples.append((cell_set_hlca_vertex, "EXPRESSES", gene_vertex))
    for gene in cellref_genes:
        gene_vertex = vertex_collections["gene"].get(gene)
        triples.append((cell_set_cellref_vertex, "EXPRESSES", gene_vertex))

    # Cell set SOURCE publication triples
    triples.extend(
        [
            (cell_set_hlca_vertex, "SOURCE", publication_hlca_vertex),
            (cell_set_cellref_vertex, "SOURCE", publication_cellref_vertex),
        ]
    )

    # Gene PART_OF biomarker combination triples
    for gene in hlca_genes:
        gene_vertex = vertex_collections["gene"].get(gene)
        triples.append((gene_vertex, "PART_OF", biomarker_combination_hlca_vertex))
    for gene in cellref_genes:
        gene_vertex = vertex_collections["gene"].get(gene)
        triples.append((gene_vertex, "PART_OF", biomarker_combination_cellref_vertex))

    for triple in triples:

        from_vertex = triple[0]
        predicate = triple[1]
        to_vertex = triple[2]

        if from_vertex == {} or to_vertex == {}:
            continue

        edge = {
            "_key": from_vertex["_key"] + "-" + to_vertex["_key"],
            "_from": from_vertex["_id"],
            "_to": to_vertex["_id"],
            "name": predicate,
        }

        from_vertex_collection = from_vertex["_id"].split("/")[0]
        to_vertex_collection = to_vertex["_id"].split("/")[0]
        edge_name = from_vertex_collection + "-" + to_vertex_collection

        if not edge_collections[edge_name].has(edge):
            edge_collections[edge_name].insert(edge)