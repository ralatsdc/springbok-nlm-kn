from datetime import datetime
from pathlib import Path
from pprint import pprint
import re
from urllib.parse import urlparse

import ArangoDB as adb
from lxml import etree
from rdflib import Graph
from rdflib.term import BNode, Literal


URIREF_PATTERN = re.compile(r"/obo/([A-Za-z]*)_([A-Z0-9]*)")
VALID_VERTICES = ["UBERON", "CL", "GO", "NCBITaxon", "PR", "PATO"]


def count_triple_types(rdf_graph):

    triple_types = {}

    for s, p, o in rdf_graph:
        triple_type = (type(s), type(p), type(o))

        if triple_type not in triple_types:
            triple_types[triple_type] = 1

        else:
            triple_types[triple_type] += 1

    return triple_types


def collect_fnode_triples(rdf_graph):

    triples = []

    for s, p, o in rdf_graph:

        if isinstance(s, BNode) or isinstance(o, BNode):
            continue

        triples.append((s, p, o))

    return triples


def parse_term(term, ols=None):

    path = urlparse(term).path
    fragment = urlparse(term).fragment

    match = URIREF_PATTERN.match(path)

    if match is not None:

        oid = match.group(1)
        if oid == "GOREL":
            print(f"Invalid Ontology ID: 'GOREL' for term: {term}")
            return None, None, None, None, None

        number = match.group(2)
        if len(oid) == 0 or len(number) == 0:
            print(f"Did not match ontology id or number for term: {term}")
            return None, None, None, None, None

        term = f"{oid}_{number}"

        if ols is not None and term in ols:
            return oid, number, term, ols[term], "class"

        else:
            return oid, number, term, None, "class"

    elif fragment != "":
        return None, None, None, fragment, "predicate"

    else:
        return None, None, None, Path(path).stem, "literal"


def collect_bnode_triple_sets(rdf_graph, triple_sets, use="subject", ols=None):

    for s, p, o in rdf_graph:

        if isinstance(s, BNode) and isinstance(o, BNode):
            continue

        if use == "subject":
            n = s

        elif use == "object":
            n = o

        else:
            raise Exception("Must use 'subject' or 'object'")

        if not isinstance(n, BNode):
            continue

        if n not in triple_sets:
            triple_sets[n] = {}
            triple_sets[n]["relation"] = []
            triple_sets[n]["annotation"] = []
            triple_sets[n]["literal"] = []
            triple_sets[n]["class"] = []
            triple_sets[n]["other"] = []

        _, _, _, _, s_term_type = parse_term(s, ols)
        _, _, _, p_fragment, _ = parse_term(p, ols)
        _, _, _, _, o_term_type = parse_term(o, ols)

        if p_fragment in ["someValuesFrom", "onProperty", "subClassOf"]:
            triple_sets[n]["relation"].append((s, p, o))

        elif p_fragment in ["annotatedSource", "annotatedProperty", "annotatedTarget"]:
            triple_sets[n]["annotation"].append((s, p, o))

        elif p_fragment in ["hasDbXref", "source"]:
            triple_sets[n]["literal"].append((s, p, o))

        elif s_term_type == "class" or o_term_type == "class":
            triple_sets[n]["class"].append((s, p, o))

        else:
            triple_sets[n]["other"].append((s, p, o))


def get_fnode(s, o):

    if isinstance(s, BNode) and isinstance(o, BNode):
        raise Exception("Both s and o are blank")

    if not isinstance(s, BNode) and not isinstance(o, BNode):
        raise Exception("Both s and o are filled")

    if isinstance(s, BNode):
        return o

    else:
        return s


def create_bnode_triples_from_bnode_triple_set(triple_set, set_type, ols=None):

    bnode_triples = []
    ignored_triples = []

    if set_type == "relation":
        s_p_fragment = "someValuesFrom"
        p_p_fragment = "onProperty"
        o_p_fragment = "subClassOf"

    elif set_type == "annotation":
        s_p_fragment = "annotatedSource"
        p_p_fragment = "annotatedProperty"
        o_p_fragment = "annotatedTarget"

    else:
        raise Exception("Set type must be 'relation' or 'annotation'")

    if len(triple_set[set_type]) == 3:

        created_s = None
        created_p = None
        created_o = None

        for s, p, o in triple_set[set_type]:

            _, _, _, p_fragment, _ = parse_term(p, ols=ols)

            if p_fragment == s_p_fragment:
                created_s = get_fnode(s, o)

            if p_fragment == p_p_fragment:
                created_p = get_fnode(s, o)

            if p_fragment == o_p_fragment:
                created_o = get_fnode(s, o)

        if created_s is not None and created_p is not None and created_o is not None:
            bnode_triples.append((created_s, created_p, created_o))

            if set_type == "annotation":

                for s, p, o in triple_set["literal"]:
                    bnode_triples.append((created_s, p, o))

        else:
            pprint(f"Invalid triple_set['{set_type}']: {triple_set[set_type]}")
            ignored_triples.extend(triple_set[set_type])

            if set_type == "annotation":
                ignored_triples.extend(triple_set["literal"])

    elif len(triple_set[set_type]) != 0:
        ignored_triples.extend(triple_set[set_type])

        if set_type == "annotation":
            ignored_triples.extend(triple_set["literal"])

    return bnode_triples, ignored_triples


def create_bnode_triples_from_bnode_triple_sets(triple_sets, ols=None):

    bnode_triples = []
    ignored_triples = []

    for bnode, triple_set in triple_sets.items():

        relation_bnode_triples, relation_ignored_triples = (
            create_bnode_triples_from_bnode_triple_set(triple_set, "relation", ols=ols)
        )

        bnode_triples.extend(relation_bnode_triples)
        ignored_triples.extend(relation_ignored_triples)

        annotation_bnode_triples, annotation_ignored_triples = (
            create_bnode_triples_from_bnode_triple_set(
                triple_set, "annotation", ols=ols
            )
        )

        bnode_triples.extend(annotation_bnode_triples)
        ignored_triples.extend(annotation_ignored_triples)

        ignored_triples.extend(triple_set["class"])
        ignored_triples.extend(triple_set["other"])

    return bnode_triples, ignored_triples


def parse_ols(ols_dir, ols_fnm):

    ols = {}
    ids = set()

    root = etree.parse(Path(ols_dir) / ols_fnm)

    owl = "{http://www.w3.org/2002/07/owl#}"
    rdf = "{http://www.w3.org/1999/02/22-rdf-syntax-ns#}"
    rdfs = "{http://www.w3.org/2000/01/rdf-schema#}"

    about_element_types = [
        "AnnotationProperty",
        "ObjectProperty",
        "DatatypeProperty",
        "Class",
        "Description",
    ]

    for about_element_type in about_element_types:

        for about_element in root.iter(f"{owl}{about_element_type}"):

            uriref = about_element.get(f"{rdf}about")
            if uriref is None:
                continue

            id, number, term, _, _ = parse_term(uriref)
            if id is None:
                continue

            label_element = about_element.find(f"{rdfs}label")
            if label_element is None:
                continue

            label = label_element.text
            ols[term] = label
            ids.add(id)

    return ols, ids


def create_or_get_vertices_from_triple(
    adb_graph, vertex_collections, s, p, o, ols=None
):

    if isinstance(o, Literal):
        return

    vertices = []
    for term in [s, o]:
        oid, number, term, fragment, term_type = parse_term(term, ols=ols)

        if term_type != "class":
            continue

        vertex_name = oid
        vertex_key = number
        vertex_term = term

        if vertex_name not in VALID_VERTICES:
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
            vertex_collections[vertex_name].insert(vertex)

        vertices.append(vertex_collections[vertex_name].get(vertex_key))

    return vertex_collections


def create_or_get_edge_from_triple(
    adb_graph, vertex_collections, edge_collections, s, p, o, ols=None
):

    if isinstance(o, Literal):
        return

    s_oid, s_number, s_term, s_fragment, s_term_type = parse_term(s, ols=ols)

    if s_term_type != "class":
        return

    from_vertex_name = s_oid
    from_vertex_key = s_number
    from_vertex_term = s_term

    if from_vertex_name not in VALID_VERTICES:
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
        vertex_collections[from_vertex_name].insert(vertex)

    p_oid, p_number, p_term, p_fragment, p_term_type = parse_term(p, ols=ols)

    if not (p_term_type == "predicate" or (p_term_type == "class" and p_fragment is not None)):
        return

    predicate = p_fragment

    o_oid, o_number, o_term, o_fragment, o_term_type = parse_term(o, ols=ols)

    if o_term_type != "class" and o_term_type != "literal":
        return

    to_vertex_name = o_oid
    to_vertex_key = o_number
    to_vertex_term = o_term

    if to_vertex_name not in VALID_VERTICES:
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
        edge_collections[edge_name].insert(edge)

    return vertex_collections, edge_collections


def update_vertex_from_triple(vertex_collections, s, p, o, ols=None):

    if not isinstance(o, Literal):
        return

    s_oid, s_number, s_term, s_fragment, s_term_type = parse_term(s, ols=ols)

    if s_term_type != "class":
        return

    vertex_name = s_oid
    vertex_key = s_number
    vertex_term = s_term

    if vertex_name not in VALID_VERTICES:
        return

    p_oid, p_number, p_term, p_fragment, p_term_type = parse_term(p, ols=ols)

    if p_term_type != "predicate":
        return

    predicate = p_fragment

    if isinstance(o.value, datetime):
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


def load_triples_into_adb_graph(
    triples, adb_graph, vertex_collections, edge_collections, ols=None
):

    for s, p, o in triples:

        create_or_get_vertices_from_triple(
            adb_graph, vertex_collections, s, p, o, ols=ols
        )
        create_or_get_edge_from_triple(
            adb_graph, vertex_collections, edge_collections, s, p, o, ols=ols
        )

    for s, p, o in triples:

        update_vertex_from_triple(vertex_collections, s, p, o, ols=ols)

    return vertex_collections, edge_collections


if __name__ == "__main__":

    ols_dir = "/Users/raymondleclair/Projects/NLM/NLM-KB/springbok-ncbi-cell/ncbi-cell/data/ols"
    ols_fnm = "ro.owl"

    ols, _ = parse_ols(ols_dir, ols_fnm)

    bioportal_dir = "/Users/raymondleclair/Projects/NLM/NLM-KB/springbok-ncbi-cell/ncbi-cell/data/bioportal"
    cl_fnm = "general_cell_types_upper_slim.owl"

    _, ids = parse_ols(bioportal_dir, cl_fnm)
    print(ids)

    rdf_graph = Graph()
    rdf_graph.parse(Path(bioportal_dir) / cl_fnm)

    triple_types = count_triple_types(rdf_graph)
    pprint(triple_types)

    fnode_triples = collect_fnode_triples(rdf_graph)

    with open("fnode_triples.txt", "w") as fp:
        for fnode_triple in fnode_triples:
            fp.write(str(fnode_triple) + "\n")

    bnode_triple_sets = {}

    collect_bnode_triple_sets(rdf_graph, bnode_triple_sets, use="subject", ols=ols)
    collect_bnode_triple_sets(rdf_graph, bnode_triple_sets, use="object", ols=ols)

    bnode_triples, ignored_triples = create_bnode_triples_from_bnode_triple_sets(
        bnode_triple_sets, ols=ols
    )

    with open("bnode_triples.txt", "w") as fp:
        for bnode_triple in bnode_triples:
            fp.write(str(bnode_triple) + "\n")

    fnode_triples.extend(bnode_triples)

    db_name = "BioPortal"
    graph_name = "CL"

    adb.delete_database(db_name)
    db = adb.create_or_get_database(db_name)

    adb.delete_graph(db, graph_name)
    adb_graph = adb.create_or_get_graph(db, graph_name)

    vertex_collections = {}
    edge_collections = {}

    load_triples_into_adb_graph(
        fnode_triples, adb_graph, vertex_collections, edge_collections, ols=ols
    )
