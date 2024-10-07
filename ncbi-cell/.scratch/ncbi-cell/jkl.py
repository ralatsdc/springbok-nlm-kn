import datetime
from pathlib import Path
from pprint import pprint
import re
from warnings import warn
from urllib.parse import urlparse

import ArangoDB as adb
from rdflib import Graph
from rdflib.term import BNode, Literal, URIRef

URIREF_PATTERN = re.compile(r"/obo/([A-Za-z]*)_([A-Z0-9]*)")
PREDICATE_CLASSES = {
    "IAO_0000028": "symbol",
    "IAO_0000115": "definition",
    "IAO_0000116": "editor note",
    "IAO_0000231": "has obsolescence reason",
    "IAO_0000233": "term tracker item",
    "IAO_0000424": "expand expression to",
    "IAO_0000700": "has ontology root term",
    "IAO_0100001": "term replaced by",
    "OMO_0002000": "defined by construct",
    "RO_0002161": "never in taxon",
    "RO_0002175": "present in taxon",
    "RO_0004050": "is negative form of",
}


def count_triple_types(rdf_graph):
    triple_types = {}
    for s, p, o in rdf_graph:
        triple_type = (type(s), type(p), type(o))
        if triple_type not in triple_types:
            triple_types[triple_type] = 1
        else:
            triple_types[triple_type] += 1
    return triple_types


def get_triples_by_type(rdf_graph, triple_type):
    triples = []
    for s, p, o in rdf_graph:
        if (type(s), type(p), type(o)) == triple_type:
            triples.append((s, p, o))
    return triples


def count_triple_components(rdf_graph):

    subjects = {}
    predicates = {}
    objects = {}

    for s, p, o in rdf_graph:

        if s not in subjects:
            subjects[s] = 1
        else:
            subjects[s] += 1

        if p not in predicates:
            predicates[p] = 1
        else:
            predicates[p] += 1

        if o not in objects:
            objects[o] = 1
        else:
            objects[o] += 1

    return subjects, predicates, objects


def parse_term(term):

    path = urlparse(term).path
    fragment = urlparse(term).fragment

    match = URIREF_PATTERN.match(path)

    if match is not None:

        oid = match.group(1)
        if oid == "GOREL":
            warn("Invalid Ontology ID: 'GOREL'")
            return tuple()
        number = match.group(2)
        if len(oid) == 0 or len(number) == 0:
            warn(f"Did not match ontology id or number for: {term}")
            return tuple()
        term = f"{oid}_{number}"

        if term in PREDICATE_CLASSES:
            return (oid, number, term, PREDICATE_CLASSES[term])

        else:
            return (oid, number, term)

    elif fragment != "":
        return (fragment,)

    else:
        return (Path(path).stem,)


def create_or_get_vertices_from_triple(adb_graph, vertex_collections, s, p, o):

    if not isinstance(s, URIRef) or isinstance(o, Literal):
        return

    # (rdflib.term.URIRef, rdflib.term.URIRef, rdflib.term.URIRef)
    # (rdflib.term.URIRef, rdflib.term.URIRef, rdflib.term.BNode)

    vertices = []
    for term in [s, o]:
        term_tuple = parse_term(term)

        if len(term_tuple) == 3:
            vertex_name, vertex_key, vertex_term = term_tuple

        elif len(term_tuple) == 1 and isinstance(term, BNode):
            vertex_name = "BNode"
            vertex_key = term_tuple[0]
            vertex_term = "none"

        else:
            continue

        if vertex_name not in vertex_collections:
            vertex_collections[vertex_name] = adb.create_or_get_vertex_collection(
                adb_graph, vertex_name
            )

        if not vertex_collections[vertex_name].has(vertex_key):
            vertex = {
                "_key": vertex_key,
                "term": vertex_term,
            }
            print(vertex)
            vertex_collections[vertex_name].insert(vertex)

        vertices.append(vertex_collections[vertex_name].get(vertex_key))

    return vertex_collections


def create_or_get_edge_from_triple(
    adb_graph, vertex_collections, edge_collections, s, p, o
):

    if isinstance(o, Literal):
        return

    # (rdflib.term.BNode, rdflib.term.URIRef, rdflib.term.URIRef): 10277,
    # (rdflib.term.URIRef, rdflib.term.URIRef, rdflib.term.URIRef): 3957,
    # (rdflib.term.BNode, rdflib.term.URIRef, rdflib.term.BNode): 1298,
    # (rdflib.term.URIRef, rdflib.term.URIRef, rdflib.term.BNode): 1166}

    s_tuple = parse_term(s)
    if len(s_tuple) == 3:
        from_vertex_name, from_vertex_key, from_vertex_term = s_tuple
    elif len(s_tuple) == 1 and isinstance(s, BNode):
        from_vertex_name = "BNode"
        from_vertex_key = s_tuple[0]
        from_vertex_term = "none"
    else:
        return
    if from_vertex_name not in vertex_collections:
        vertex_collections[from_vertex_name] = adb.create_or_get_vertex_collection(
            adb_graph, from_vertex_name
        )
    if not vertex_collections[from_vertex_name].has(from_vertex_key):
        vertex = {
            "_key": from_vertex_key,
            "term": from_vertex_term,
        }
        print(vertex)
        vertex_collections[from_vertex_name].insert(vertex)

    p_tuple = parse_term(p)
    if len(p_tuple) == 1:
        predicate = p_tuple[0]
    elif len(p_tuple) == 4:
        predicate = p_tuple[3]
    else:
        return

    o_tuple = parse_term(o)
    if len(o_tuple) == 3:
        to_vertex_name, to_vertex_key, to_vertex_term = o_tuple
    elif len(o_tuple) == 1 and isinstance(o, BNode):
        to_vertex_name = "BNode"
        to_vertex_key = o_tuple[0]
        to_vertex_term = "none"
    else:
        return
    if to_vertex_name not in vertex_collections:
        vertex_collections[to_vertex_name] = adb.create_or_get_vertex_collection(
            adb_graph, to_vertex_name
        )
    if not vertex_collections[to_vertex_name].has(to_vertex_key):
        vertex = {
            "_key": to_vertex_key,
            "term": to_vertex_term,
        }
        print(vertex)
        vertex_collections[to_vertex_name].insert(vertex)

    edge_name = f"{from_vertex_name}-{to_vertex_name}"
    edge_key = f"{from_vertex_key}-{to_vertex_key}"

    if edge_name not in edge_collections:
        edge_collections[edge_name] = adb.create_or_get_edge_collection(
            adb_graph, from_vertex_name, to_vertex_name
        )[0]

    if not edge_collections[edge_name].has(edge_key):
        edge = {
            "_key": edge_key,
            "_from": f"{from_vertex_name}/{from_vertex_key}",
            "_to": f"{to_vertex_name}/{to_vertex_key}",
            "label": predicate,
        }
        print(edge)
        edge_collections[edge_name].insert(edge)

    return vertex_collections, edge_collections


def update_vertex_from_triple(vertex_collections, s, p, o):

    if not isinstance(s, URIRef) or not isinstance(o, Literal):
        return

    # (rdflib.term.URIRef, rdflib.term.URIRef, rdflib.term.URIRef)
    # (rdflib.term.URIRef, rdflib.term.URIRef, rdflib.term.BNode)

    s_tuple = parse_term(s)
    p_tuple = parse_term(p)

    if not (len(s_tuple) == 3 and (len(p_tuple) == 1 or len(p_tuple) == 4)):
        return

    vertex_name, vertex_key, _ = s_tuple

    if len(p_tuple) == 1:
        predicate = p_tuple[0]
    elif len(p_tuple) == 4:
        predicate = p_tuple[3]
    else:
        return

    if isinstance(o.value, datetime.datetime):
        value = str(o.value)
    else:
        value = o.value

    if vertex_name in vertex_collections and vertex_collections[vertex_name].has(
        vertex_key
    ):
        vertex = vertex_collections[vertex_name].get(vertex_key)
        if not predicate in vertex:
            vertex[predicate] = value
        else:
            if not isinstance(vertex[predicate], list):
                vertex[predicate] = [vertex[predicate]]
            vertex[predicate].append(value)
        vertex_collections[vertex_name].update(vertex)
        return vertex

    else:
        return


def load_rdf_graph_into_adb_graph(
    rdf_graph, adb_graph, vertex_collections, edge_collections, term=None
):

    for s, p, o in rdf_graph:

        if term is not None and s != term:
            continue

        create_or_get_vertices_from_triple(adb_graph, vertex_collections, s, p, o)
        create_or_get_edge_from_triple(
            adb_graph, vertex_collections, edge_collections, s, p, o
        )

    for s, p, o in rdf_graph:

        if term is not None and s != term:
            continue

        update_vertex_from_triple(vertex_collections, s, p, o)

    return vertex_collections, edge_collections


def print_triples_with_term(rdf_graph, term):

    triple_types = count_triple_types(rdf_graph)

    print()
    pprint(triple_types)

    # subjects, predicates, objects = count_triple_components(rdf_graph)

    bnodes = []
    for s, p, o in rdf_graph:
        if s == term and isinstance(o, BNode):
            print()
            print(f"s: {s}")
            print(f"p: {p}")
            print(f"o: {o}")
            bnodes.append(o)

    print()
    print("bnodes")
    for s, p, o in rdf_graph:
        if s in bnodes:
            print()
            print(f"s: {s}")
            print(f"p: {p}")
            print(f"o: {o}")


if __name__ == "__main__":

    bioportal_dir = Path(
        "/Users/raymondleclair/Projects/NLM/NLM-KB/springbok-ncbi-cell/ncbi-cell/data/bioportal"
    )
    # cl_fnm = "cl.owl"
    cl_fnm = "general_cell_types_upper_slim.owl"

    rdf_graph = Graph()
    rdf_graph.parse(str(bioportal_dir / cl_fnm))

    adb.delete_database("BioPortal")
    db = adb.create_or_get_database("BioPortal")
    adb_graph = adb.create_or_get_graph(db, "CL")

    vertex_collections = {}
    edge_collections = {}

    # print_triples_with_term(rdf_graph, URIRef("http://purl.obolibrary.org/obo/CL_0000235"))

    load_rdf_graph_into_adb_graph(
        rdf_graph,
        adb_graph,
        vertex_collections,
        edge_collections,
        term=URIRef("http://purl.obolibrary.org/obo/CL_0000235"),
    )

    # input("Hit return to continue?")

    do_add_bnodes = True
    while do_add_bnodes:
        for bnode_vertex in vertex_collections["BNode"]:
            load_rdf_graph_into_adb_graph(
                rdf_graph,
                adb_graph,
                vertex_collections,
                edge_collections,
                term=BNode(bnode_vertex["_key"]),
            )
        # do_add_bnodes = input("Continue? [Y/n]: ") != "n"
        do_add_bnodes = False


    for edge_name, edge_collection in edge_collections.items():
        if 'BNode' not in edge_name or edge_name == 'BNode-BNode':
            continue
        for edge in edge_collection:
            print(edge)

    # print_triples_with_term(rdf_graph, BNode("Nf3b88ff808404a80a0b7f29036ba8fa7"))
