from pathlib import Path
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


def count_triple_types(graph):
    triple_types = {}
    for s, p, o in graph:
        triple_type = (type(s), type(p), type(o))
        if triple_type not in triple_types:
            triple_types[triple_type] = 1
        else:
            triple_types[triple_type] += 1
    return triple_types


def get_triples_by_type(graph, triple_type):
    triples = []
    for s, p, o in graph:
        if (type(s), type(p), type(o)) == triple_type:
            triples.append((s, p, o))
    return triples


def count_triple_components(graph):

    subjects = {}
    predicates = {}
    objects = {}

    for s, p, o in graph:

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


def urirefparse(uriref):

    path = urlparse(uriref).path
    fragment = urlparse(uriref).fragment

    match = URIREF_PATTERN.match(path)

    if match is not None:

        oid = match.group(1)
        if oid == "GOREL":
            warn("Invalid Ontology ID: 'GOREL'")
            return tuple()
        number = match.group(2)
        if len(oid) == 0 or len(number) == 0:
            warn(f"Did not match ontology id or number for: {uriref}")
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


def create_or_get_vertex_from_class_tuple(graph, vertex_collections, class_tuple):

    vertex_name, vertex_key, vertex_term = class_tuple

    if vertex_name not in vertex_collections:
        vertex_collections[vertex_name] = adb.create_or_get_vertex_collection(
            graph, vertex_name
        )

    if not vertex_collections[vertex_name].has(vertex_key):
        vertex = {
            "_key": vertex_key,
            "term": vertex_term,
        }
        print(vertex)
        vertex_collections[vertex_name].insert(vertex)

    return vertex_collections[vertex_name].get(vertex_key)


def create_or_get_edge_from_class_tuples(
    graph, edge_collections, s_tuple, p_tuple, o_tuple
):

    from_vertex_name, from_vertex_key, _ = s_tuple

    predicate = ""
    if len(p_tuple) == 1:
        predicate = p_tuple[0]
    elif len(p_tuple) == 4:
        predicate = p_tuple[3]
    else:
        warn("Predicate tuple does not have one or four components")

    to_vertex_name, to_vertex_key, _ = o_tuple

    edge_name = f"{from_vertex_name}-{to_vertex_name}"
    edge_key = f"{from_vertex_key}-{to_vertex_key}"

    if edge_name not in edge_collections:
        edge_collections[edge_name] = adb.create_or_get_edge_collection(
            graph, from_vertex_name, to_vertex_name
        )[0]

    if not edge_collections[edge_name].has(edge_key):
        edge = {
            "_key": edge_key,
            "_from": f"{from_vertex_name}/{from_vertex_key}",
            "_to": f"{to_vertex_name}/{to_vertex_key}",
            "label": predicate,
        }
        print(s_tuple, p_tuple, o_tuple)
        print(edge)
        edge_collections[edge_name].insert(edge)

    return edge_collections[edge_name].get(edge_key)


if __name__ == "__main__":

    bioportal_dir = Path(
        "/Users/raymondleclair/Projects/NLM/NLM-KB/springbok-ncbi-cell/ncbi-cell/data/bioportal"
    )
    # cl_fnm = "cl.owl"
    cl_fnm = "general_cell_types_upper_slim.owl"

    rdf_graph = Graph()
    rdf_graph.parse(str(bioportal_dir / cl_fnm))

    triple_types = count_triple_types(rdf_graph)
    subjects, predicates, objects = count_triple_components(rdf_graph)

    class_triples = get_triples_by_type(rdf_graph, (URIRef, URIRef, URIRef))
    label_triples = get_triples_by_type(rdf_graph, (URIRef, URIRef, Literal))

    adb.delete_database("BioPortal")
    db = adb.create_or_get_database("BioPortal")
    adb_graph = adb.create_or_get_graph(db, "CL")

    vertex_collections = {}
    edge_collections = {}

    for s, p, o in class_triples:

        s_tuple = urirefparse(s)
        p_tuple = urirefparse(p)
        o_tuple = urirefparse(o)

        if (
            len(s_tuple) == 3
            and (len(p_tuple) == 1 or len(p_tuple) == 4)
            and len(o_tuple) == 3
        ):

            s_vertex = create_or_get_vertex_from_class_tuple(
                adb_graph, vertex_collections, s_tuple
            )

            o_vertex = create_or_get_vertex_from_class_tuple(
                adb_graph, vertex_collections, o_tuple
            )

            s_o_edge = create_or_get_edge_from_class_tuples(
                adb_graph, edge_collections, s_tuple, p_tuple, o_tuple
            )
