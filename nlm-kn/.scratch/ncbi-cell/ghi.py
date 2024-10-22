from pathlib import Path
import re
from urllib.parse import urlparse
from warnings import warn

from rdflib import Graph
from rdflib.term import BNode, Literal, URIRef

BIOPORTAL_DIR = Path(
    "/Users/raymondleclair/Projects/NLM/NLM-KB/springbok-ncbi-cell/ncbi-cell/data/bioportal"
)

# CL_FNM = "cl.owl"
CL_FNM = "general_cell_types_upper_slim.owl"

GRAPH = Graph()
GRAPH.parse(str(BIOPORTAL_DIR / CL_FNM))

URIREF_PATTERN = re.compile(r"/obo/([A-Z]*)_([A-Z0-9]*)")


def urirefparse(uriref):
    path = urlparse(uriref).path

    match = URIREF_PATTERN.match(path)
    if match is None:
        warn("URIRef pattern did not match")
        return None, None, None

    oid = match.group(1)
    number = match.group(2)
    if len(oid) == 0 or len(number) == 0:
        warn("Did not match id or number")
        return None, None, None

    term = f"{oid}_{number}"

    return oid, number, term


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


def count_spo(graph):

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


def count_predicates(graph):
    predicates = {}
    for _s, p, _o in graph:
        if p not in predicates:
            predicates[p] = 1
        else:
            predicates[p] += 1
    return predicates


def get_class_triples_for_not_predicate(graph, predicate):
    triples = []
    for s, p, o in graph:
        if p != predicate and (isinstance(s, URIRef) and isinstance(o, URIRef)):
            s_id, _s_number, _s_term = urirefparse(s)
            o_id, _o_number, _o_term = urirefparse(o)
            if s_id is not None and o_id is not None:
                triples.append((s, p, o))
    return triples

def get_class_triples_for_predicate(graph, predicate):
    triples = []
    for s, p, o in graph:
        if p == predicate and (isinstance(s, URIRef) and isinstance(o, URIRef)):
            s_id, _s_number, _s_term = urirefparse(s)
            o_id, _o_number, _o_term = urirefparse(o)
            if s_id is not None and o_id is not None:
                triples.append((s, p, o))
    return triples

def get_triples_for_precdicate(graph, predicate):
    triples = []
    for s, p, o in graph:
        if p == predicate:
            triples.append((s, p, o))
    return triples

subjects, predicates, objects = count_spo(GRAPH)

triple_types = count_triple_types(GRAPH)

# class_triples = get_class_triples_for_not_predicate(GRAPH, URIRef("http://www.w3.org/2000/01/rdf-schema#subClassOf"))
# subclass_of_triples = get_class_triples_for_predicate(GRAPH, URIRef("http://www.w3.org/2000/01/rdf-schema#subClassOf"))
# label_triples = get_triples_for_precdicate(GRAPH, URIRef("http://www.w3.org/2000/01/rdf-schema#label"))
# definition_triples = get_triples_for_precdicate(GRAPH, URIRef("http://purl.obolibrary.org/obo/IAO_0000115"))
