import argparse
from pathlib import Path

import pandas as pd

import ArangoDB as adb
from CellOntology import (
    load_triples_into_adb_graph,
    parse_ols,
    parse_term,
    VALID_VERTICES,
)


def read_excel(nlm_kn_dirname, nlm_kn_filename):
    """Read NLM KN Ontology Excel file.

    Parameters
    ----------
    nlm_kn_dirname : str | Path
        The name of directory containing NLM KN ontology
    nlm_kn_filename : str
        The name of file containing NLM KN ontology

    Returns
    -------
    schema : pd.DataFrame
        The DataFrame containing the schema name triples
    relations : pd.DataFrame
        The DataFrame containing the relations of the schema names to
        CURIEs
    """
    schema = pd.read_excel(Path(nlm_kn_dirname) / nlm_kn_filename, sheet_name=0).iloc[
        :, 0:3
    ]
    relations = pd.read_excel(
        Path(nlm_kn_dirname) / nlm_kn_filename, sheet_name=1
    ).iloc[:, 0:3]

    return schema, relations


def create_triples(schema, relations, ro=None):
    """Create triples by translating schema names using CURIEs and a default URL.

    Parameters
    ----------
    schema : pd.DataFrame
        The DataFrame containing the schema name triples
    relations : pd.DataFrame
        The DataFrame containing the relations of the schema names to
        CURIEs

    Returns
    -------
    triples : list(tuple)
        List of tuples which contain each triple
    """
    triples = []

    # Create mapping from name to CURIE
    n2c = {}
    for _, row in relations.iterrows():
        n2c[row.iloc[0]] = row.iloc[2]

    # Create triples by translating names using CURIEs and a default
    # URL
    ids = set()
    for _, row in schema.iterrows():
        triple = []
        for item in row:
            if not item in n2c:
                print(f"Skipping item: {item}")
                continue
            term = f"http://purl.obolibrary.org/obo/{n2c[item].replace(':', '_')}"
            triple.append(term)
            id, _, _, _, _ = parse_term(term, ro=ro)
            ids.add(id)
        if len(triple) != 3:
            print(f"Skipping row: {row.to_list()}")
            continue
        triples.append(tuple(triple))

    return triples, ids


def main():

    parser = argparse.ArgumentParser(description="Load NLM-KN Ontology")
    parser.add_argument(
        "--nlm-kn-dirname",
        default=Path("../data/nlm-kn"),
        help="Name of directory containing NLM KN ontology",
    )
    parser.add_argument(
        "--nlm-kn-filename",
        default=Path("Cell_Phenotype_KG_Schema_v2.xlsx"),
        help="Name of file containing NLM KN ontology",
    )
    group = parser.add_argument_group("Cell Ontology (CL)", "Version of the CL loaded")
    exclusive_group = group.add_mutually_exclusive_group(required=True)
    exclusive_group.add_argument(
        "--slim", action="store_true", help="Slim ontology loaded"
    )
    exclusive_group.add_argument(
        "--full", action="store_true", help="Full ontology loaded"
    )
    group.add_argument("--include-bnodes", action="store_true", help="BNodes included")
    parser.add_argument(
        "--ols-dirname",
        default=Path("../data/ols"),
        help="Name of directory containing ontologies downloaded from the Ontology Lookup Service",
    )

    args = parser.parse_args()

    if args.slim:
        cl_filename = "general_cell_types_upper_slim.owl"
        db_name = "BioPortal-Slim"
        graph_name = "CL-Slim"

    if args.full:
        cl_filename = "cl.owl"
        db_name = "BioPortal-Full"
        graph_name = "CL-Full"

    if args.include_bnodes:
        db_name += "-BNodes"
        graph_name += "-BNodes"

    db_name += "-Test"
    graph_name += "-Test"

    ro_filename = "ro.owl"
    log_filename = f"{graph_name}.log"

    print(f"Reading {args.nlm_kn_dirname / args.nlm_kn_filename}")
    schema, relations = read_excel(args.nlm_kn_dirname, args.nlm_kn_filename)

    print("Creating triples")
    ro, _, _ = parse_ols(args.ols_dirname, ro_filename)
    triples, ids = create_triples(schema, relations, ro=ro)

    # print("Creating ArangoDB database and graph, and loading triples")
    # VALID_VERTICES.update(ids.difference(set([None, "BFO", "IAO", "RO"])))
    # db = adb.create_or_get_database(db_name)
    # adb_graph = adb.create_or_get_graph(db, graph_name)
    # vertex_collections = {}
    # edge_collections = {}
    # load_triples_into_adb_graph(
    #     triples, adb_graph, vertex_collections, edge_collections, ro=ro
    # )


if __name__ == "__main__":
    main()
