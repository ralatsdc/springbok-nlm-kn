from pathlib import Path
import pprint
import re
from urllib.parse import urlparse
from warnings import warn

from rdflib import Graph
from rdflib.term import BNode, Literal, URIRef
from lxml import etree

import ArangoDB as adb


vertex_names = ["CARO", "CHEBI", "CL", "GO", "PATO", "PR", "UBERON"]

from_to_vertex_name_sets = [
    ("CHEBI", "BFO"),
    ("CHEBI", "PR"),
    ("CL", "BFO"),
    ("CL", "CL"),
    ("CL", "GO"),
    ("CL", "PR"),
    ("CL", "UBERON"),
    ("GO", "BFO"),
    ("GO", "CARO"),
    ("GO", "GO"),
    ("GO", "PR"),
    ("GO", "UBERON"),
    ("IAO", "BFO"),
    ("IAO", "IAO"),
    ("PATO", "BFO"),
    ("PATO", "PATO"),
    ("PR", "BFO"),
    ("PR", "CHEBI"),
    ("PR", "GO"),
    ("PR", "PR"),
    ("RO", "BFO"),
    ("UBERON", "BFO"),
    ("UBERON", "GO"),
    ("UBERON", "PR"),
    ("UBERON", "UBERON"),
]

vertices = {}
edges = {}

vertex_collections = {}
edge_collections = {}

adb.delete_database("BioPortal")
db = adb.create_or_get_database("BioPortal")
graph = adb.create_or_get_graph(db, "CL")

for vertex_name in vertex_names:
    vertices[vertex_name] = []
    vertex_collections[vertex_name] = adb.create_or_get_vertex_collection(
        graph, vertex_name
    )

for from_vertex_name, to_vertex_name in from_to_vertex_name_sets:
    edge_name = f"{from_vertex_name}-{to_vertex_name}"
    edges[edge_name] = []
    edge_collections[edge_name], _ = adb.create_or_get_edge_collection(
        graph, from_vertex_name, to_vertex_name
    )


uriref_pattern = re.compile(r"/obo/([A-Z]*)_([A-Z0-9]*)")


def urirefparse(uriref):
    path = urlparse(uriref).path

    match = uriref_pattern.match(path)
    if match is None:
        warn("URIRef pattern did not match")
        return None, None, None

    id = match.group(1)
    number = match.group(2)
    if len(id) == 0 or len(number) == 0:
        warn("Did not match id or number")
        return None, None, None

    term = f"{id}_{number}"

    return id, number, term


base = "{http://purl.obolibrary.org/obo/cl.owl}"
cl = "{http://purl.obolibrary.org/obo/cl#}"
dc = "{http://purl.org/dc/elements/1.1/}"
go = "{http://purl.obolibrary.org/obo/go#}"
pr = "{http://purl.obolibrary.org/obo/pr#}"
obo = "{http://purl.obolibrary.org/obo/}"
owl = "{http://www.w3.org/2002/07/owl#}"
rdf = "{http://www.w3.org/1999/02/22-rdf-syntax-ns#}"
xml = "{http://www.w3.org/XML/1998/namespace}"
xsd = "{http://www.w3.org/2001/XMLSchema#}"
core = "{http://purl.obolibrary.org/obo/uberon/core#}"
foaf = "{http://xmlns.com/foaf/0.1/}"
pato = "{http://purl.obolibrary.org/obo/pato#}"
rdfs = "{http://www.w3.org/2000/01/rdf-schema#}"
sssom = "{https://w3id.org/sssom/}"
terms = "{http://purl.org/dc/terms/}"
uberon = "{http://purl.obolibrary.org/obo/uberon#}"
ubprop = "{http://purl.obolibrary.org/obo/ubprop#}"
subsets = "{http://purl.obolibrary.org/obo/ro/subsets#}"
oboInOwl = "{http://www.geneontology.org/formats/oboInOwl#}"
ncbitaxon = "{http://purl.obolibrary.org/obo/ncbitaxon#}"


bioportal_dir = Path(
    "/Users/raymondleclair/Projects/NLM/NLM-KB/springbok-ncbi-cell/ncbi-cell/data/bioportal"
)

cl_full_fnm = "cl.owl"
cl_slim_fnm = "general_cell_types_upper_slim.owl"

root = etree.parse(bioportal_dir / cl_slim_fnm)


ids = set()

for class_elm in root.iter(f"{owl}Class"):

    uriref = class_elm.get(f"{rdf}about")
    if uriref is None:
        continue

    id, number, term = urirefparse(uriref)
    if id is None:
        continue

    ids.add(id)
    if id in vertex_names:

        label = ""
        label_elm = class_elm.find(f"{rdfs}label")
        if label_elm is not None:

            label = label_elm.text

        d = {"_key": number, "term": term, "uriref": uriref, "label": label}

        vertices[id].append(d)

        # vertex_collections[id].insert(d)


graph = Graph()
graph.parse(str(bioportal_dir / cl_slim_fnm))

type_s = set()
type_p = set()
type_o = set()

bnode_s = []
bnode_o = []

bnode_t = {}

unique_p = set()

count_p = {}

edge_keys = set()

for s, p, o in graph:
    type_s.add(type(s))
    type_p.add(type(p))
    type_o.add(type(o))

    unique_p.add(p)

    if p not in count_p:
        count_p[p] = 1
    else:
        count_p[p] += 1

    if isinstance(s, BNode):
        bnode_s.append(s)

        if s not in bnode_t:
            bnode_t[s] = []
        else:
            bnode_t[s].append((s, p, o))

    if isinstance(o, BNode):
        bnode_o.append(o)

    if (
        isinstance(s, URIRef)
        and p == URIRef("http://www.w3.org/2000/01/rdf-schema#subClassOf")
        and isinstance(o, URIRef)
    ):

        s_id, s_number, s_term = urirefparse(s)
        o_id, o_number, o_term = urirefparse(o)
        if s_id is None or o_id is None:
            continue

        edge_keys.add((s_id, o_id))

        d = {
            "_key": f"{s_term}-{o_term}",
            "_from": f"{s_id}/{s_number}",
            "_to": f"{o_id}/{o_number}",
            "predicate": "IS_A",
        }

        edges[f"{s_id}-{o_id}"].append(d)

        # edge_collections[f"{s_id}-{o_id}"].insert(d)

# print(type_s)
# print(type_p)
# print(type_o)

# for p in unique_p:
#     print(p)

# pprint.pprint(count_p)
