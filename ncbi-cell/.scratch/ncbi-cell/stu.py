from pathlib import Path
from pprint import pprint

from rdflib import Graph
from rdflib.term import BNode, Literal, URIRef

from jkl import parse_term, PREDICATE_CLASSES


def collect_fnode_triples(graph):

    triples = []

    for s, p, o in graph:

        if isinstance(s, BNode) or isinstance(o, BNode):
            continue

        triples.append((s, p, o))

    return triples


def collect_bnode_triple_sets(graph, triple_sets, use="subject"):

    for s, p, o in graph:

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

        _, _, _, _, s_term_type = parse_term(s)
        _, _, _, p_fragment, _ = parse_term(p)
        _, _, _, _, o_term_type = parse_term(o)

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


def create_bnode_triples_from_bnode_triple_set(triple_set, set_type):

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

    bnode_triples = []
    ignored_triples = []

    if len(triple_set[set_type]) == 3:

        created_s = None
        created_p = None
        created_o = None

        for s, p, o in triple_set[set_type]:

            _, _, _, p_fragment, _ = parse_term(p)

            if p_fragment == s_p_fragment:
                created_s = get_fnode(s, o)

            if p_fragment == p_p_fragment:
                created_p = get_fnode(s, o)

            if p_fragment == o_p_fragment:
                created_o = get_fnode(s, o)

        if created_s and created_p and created_o:
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


def create_bnode_triples_from_bnode_triple_sets(triple_sets):

    bnode_triples = []
    ignored_triples = []

    for bnode, triple_set in triple_sets.items():

        relation_bnode_triples, relation_ignored_triples = (
            create_bnode_triples_from_bnode_triple_set(triple_set, "relation")
        )

        bnode_triples.extend(relation_bnode_triples)
        ignored_triples.extend(relation_ignored_triples)

        annotation_bnode_triples, annotation_ignored_triples = (
            create_bnode_triples_from_bnode_triple_set(triple_set, "annotation")
        )

        bnode_triples.extend(annotation_bnode_triples)
        ignored_triples.extend(annotation_ignored_triples)

        ignored_triples.extend(triple_set["class"])
        ignored_triples.extend(triple_set["other"])

    return bnode_triples, ignored_triples


if __name__ == "__main__":

    bioportal_dir = Path(
        "/Users/raymondleclair/Projects/NLM/NLM-KB/springbok-ncbi-cell/ncbi-cell/data/bioportal"
    )

    ontologies = {
        "CL-OWL": {
            "filename": "general_cell_types_upper_slim.owl",
            "dbname": "BioPortal-OWL",
            "graphname": "CL-OWL",
        },  # Parsed
        # "CL-RDF": {
        #     "filename": "general_cell_types_upper_slim.rdf",
        #     "dbname": "BioPortal-RDF",
        #     "graphname": "CL-RDF",
        # },  # Parsed
        # "CL-INFERRED-OWL": {
        #     "filename": "general_cell_types_upper_slim_inferred.owl",
        #     "dbname": "BioPortal-INFERRED-OWL",
        #     "graphname": "CL-INFERRED-OWL",
        # },  # Did not parse
        # "CL-INFERRED-RDF": {
        #     "filename": "general_cell_types_upper_slim_inferred.rdf",
        #     "dbname": "BioPortal-INFERRED-RDF",
        #     "graphname": "CL-INFERRED-RDF",
        # },  # Parsed
    }

    for k, v in ontologies.items():
        v["graph"] = Graph()
        try:
            v["graph"].parse(str(bioportal_dir / v["filename"]))
        except Exception as exc:
            print(f"Could not parse {v['filename']}: {exc}")

        fnode_triples = collect_fnode_triples(v["graph"])

        bnode_triple_sets = {}

        collect_bnode_triple_sets(v["graph"], bnode_triple_sets, use="subject")
        collect_bnode_triple_sets(v["graph"], bnode_triple_sets, use="object")

        bnode_triples, ignored_triples = create_bnode_triples_from_bnode_triple_sets(
            bnode_triple_sets
        )
