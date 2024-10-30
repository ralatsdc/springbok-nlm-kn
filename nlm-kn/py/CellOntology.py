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
VALID_VERTICES = ["UBERON", "CL", "GO", "NCBITaxon", "PR", "PATO", "CHEBI"]


def parse_ols(ols_dir, ols_fnm):
    """(1)Parse ontology XML downloaded from the Ontology Lookup Service
    to create a mapping from term to label, from label to term, and a
    list of ontology identifiers.

    Parameters
    ----------
    ols_dir : str
        Name of directory containing downloaded ontology XML
    ols_fnm : str
        Name of downloaded ontology XML file

    Returns
    -------
    t2l : dict
        Dictionary mapping ontology term to label
    l2t : dict
        Dictionary mapping ontology label to term
    ids : set
        Set containing all ontology identifiers found
    """
    root = etree.parse(Path(ols_dir) / ols_fnm)

    # Define OWL XML namespaces and element types expected to contain
    # an RDF about attribute
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

    t2l = {}
    ids = set()

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

            # Collect arsed an ontology identifier, term, and label
            t2l[term] = label
            ids.add(id)

    # Invert the term to label dictionary
    l2t = {v: k for k, v in t2l.items()}

    return t2l, l2t, ids


def parse_term(term, ro=None):
    """(2)Parse an rdflib term first as an URIRef that identifies a
    class, including relationship classes, then a predicate, BNode, or
    Literal.

    Parameters
    ----------
    term : rdflib.term.* | str
        An rdflib term: BNode, Literal, or URIRef or equivalent string
    ro : None | dict
        A dictionary mapping relationship ontology terms to labels

    Returns
    -------
    tuple
        Contains ontology identifier, number, and term, label or
        literal value, and type ('class', 'predicate', or 'literal'),
        in which any element may also be None
    """
    path = urlparse(term).path
    fragment = urlparse(term).fragment

    match = URIREF_PATTERN.match(path)

    if match is not None:

        oid = match.group(1)
        if oid == "GOREL":
            # Identifier not found in the Ontology Lookup Service
            print(f"Invalid Ontology ID: 'GOREL' for term: {term}")
            return None, None, None, None, None

        number = match.group(2)
        if len(oid) == 0 or len(number) == 0:
            print(f"Did not match ontology id or number for term: {term}")
            return None, None, None, None, None

        term = f"{oid}_{number}"

        if ro is not None and term in ro:
            # Lookup label for relationship ontology term
            return oid, number, term, ro[term], "class"

        else:
            return oid, number, term, None, "class"

    elif fragment != "":
        return None, None, None, fragment, "predicate"

    elif isinstance(term, BNode):
        # Create pseudo ontology identifier, number, and term
        oid = "BNode"
        number = Path(path).stem
        term = f"{oid}_{number}"
        return oid, number, term, None, "class"

    else:
        return None, None, None, Path(path).stem, "literal"


def count_triple_types(rdf_graph):
    """(3)Count rdflib triple types, triples containing BNode,
    Literal, and URIRef.

    Parameters
    ----------
    rdf_graph : rdflib.graph.Graph
        Graph parsed by rdflib

    Returns
    -------
    triple_types : dict
        Dictionary of counts by triple type
    """
    triple_types = {}

    for s, p, o in rdf_graph:
        triple_type = (type(s), type(p), type(o))

        if triple_type not in triple_types:
            triple_types[triple_type] = 1

        else:
            triple_types[triple_type] += 1

    return triple_types


def collect_fnode_triples(rdf_graph):
    """(4)Collect filled node triples, that is, triples in which
    neither the subject nor object is a BNode, from an rdflib graph.

    Parameters
    ----------
    rdf_graph : rdflib.graph.Graph
        Graph parsed by rdflib

    Returns
    -------
    triples : list(tuples)
        List of tuples which contain each triple
    """
    triples = []

    for s, p, o in rdf_graph:

        if isinstance(s, BNode) or isinstance(o, BNode):
            continue

        triples.append((s, p, o))

    return triples


def collect_bnode_triple_sets(rdf_graph, triple_sets, use="subject", ro=None):
    """(5)Collect sets of triples each sharing a common BNode. Sets
    contain triples relating to a relation between classes, an
    annotation of a class, or a yet to be understood
    purpose. Predicate fragments, then subject or object type, are
    used to identify set type.

    Parameters
    ----------
    rdf_graph : rdflib.graph.Graph
        Graph parsed by rdflib
    triple_sets : dict
        Dictionary containing sets of triples each sharing a common
        BNode. Sets contain triples relating to a relation between
        classes, an annotation of a class, or a yet to be understood
        purpose
    use : str
        Part of triple to use as key in the dictionary ('subject', or
        'object')
    ro : None | dict
        A dictionary mapping relationship ontology terms to labels

    Returns
    -------
    None
    """
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

        _, _, _, _, s_term_type = parse_term(s, ro=ro)
        _, _, _, p_fragment, _ = parse_term(p, ro=ro)
        _, _, _, _, o_term_type = parse_term(o, ro=ro)

        # First use predicate fragments, then subject or object type,
        # to identify set type
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


def create_bnode_triples_from_bnode_triple_sets(triple_sets, ro=None):
    """(6)Create 'relation' and 'annotation' triples from BNode triple
    sets, collecting all ignored triple.

    Parameters
    ----------
    triple_sets : dict
        Dictionary containing sets of triples each sharing a common
        BNode. Sets contain triples relating to a relation between
        classes, an annotation of a class, or a yet to be understood
        purpose
    ro : None | dict
        A dictionary mapping relationship ontology terms to labels

    Returns
    -------
    bnode_triples : list(tuple)
        List of tuples which contain all 'relation' and 'annotation'
        triples
    ignored_triples : list(tuple)
        List of tuples which contain all ignored triples
    None
    """
    bnode_triples = []
    ignored_triples = []
    for bnode, triple_set in triple_sets.items():

        # Create and collect 'relation' triples
        relation_bnode_triples, relation_ignored_triples = (
            create_bnode_triples_from_bnode_triple_set(triple_set, "relation", ro=ro)
        )
        bnode_triples.extend(relation_bnode_triples)

        # Create and collect 'annotation' triples
        annotation_bnode_triples, annotation_ignored_triples = (
            create_bnode_triples_from_bnode_triple_set(triple_set, "annotation", ro=ro)
        )
        bnode_triples.extend(annotation_bnode_triples)

        # Collect ignored triples
        ignored_triples.extend(relation_ignored_triples)
        ignored_triples.extend(annotation_ignored_triples)
        ignored_triples.extend(triple_set["class"])
        ignored_triples.extend(triple_set["other"])

    return bnode_triples, ignored_triples


def create_bnode_triples_from_bnode_triple_set(triple_set, set_type, ro=None):
    """(7)Create tiples from a triple set by identifying, based on set
    type, the triple in the set that defines the subject, predicate,
    and object.

    Parameters
    ----------
    triple_sets : dict
        Dictionary containing sets of triples each sharing a common
        BNode. Sets contain triples relating to a relation between
        classes, an annotation of a class, or a yet to be understood
        purpose
    set_type : str
        Set type, ('relation' or 'annotation')
    ro : None | dict
        A dictionary mapping relationship ontology terms to labels

    Returns
    -------
    bnode_triples : list(tuple)
        List of tuples which contain all created triples
    ignored_triples : list(tuple)
        List of tuples which contain all ignored triples
    """
    bnode_triples = []
    ignored_triples = []

    # Define the fragments which identify the triple which defines the
    # subject, predicate, and object by set type
    if set_type == "relation":
        s_p_fragment = "subClassOf"
        p_p_fragment = "onProperty"
        o_p_fragment = "someValuesFrom"

    elif set_type == "annotation":
        s_p_fragment = "annotatedSource"
        p_p_fragment = "annotatedProperty"
        o_p_fragment = "annotatedTarget"

    else:
        raise Exception("Set type must be 'relation' or 'annotation'")

    # Expect exactly three triples in a set to create a triple
    if len(triple_set[set_type]) == 3:

        # Attempt to create the subject, predicate, and object
        created_s = None
        created_p = None
        created_o = None
        for s, p, o in triple_set[set_type]:

            _, _, _, p_fragment, _ = parse_term(p, ro=ro)

            if p_fragment == s_p_fragment:
                created_s = get_fnode(s, o)

            if p_fragment == p_p_fragment:
                created_p = get_fnode(s, o)

            if p_fragment == o_p_fragment:
                created_o = get_fnode(s, o)

        if created_s is not None and created_p is not None and created_o is not None:

            # Created a valid triple, so append it
            bnode_triples.append((created_s, created_p, created_o))

            if set_type == "annotation":

                # Annotation triple sets identify a class to which
                # 'literal' triple sets provide additional annotation
                for s, p, o in triple_set["literal"]:
                    bnode_triples.append((created_s, p, o))

        else:

            # Collect all invalid triple sets
            pprint(f"Invalid triple_set['{set_type}']: {triple_set[set_type]}")
            ignored_triples.extend(triple_set[set_type])

            if set_type == "annotation":

                ignored_triples.extend(triple_set["literal"])

    elif len(triple_set[set_type]) != 0:

        # Collect all invalid triple sets
        pprint(f"Invalid triple_set['{set_type}']: {triple_set[set_type]}")
        ignored_triples.extend(triple_set[set_type])

        if set_type == "annotation":
            ignored_triples.extend(triple_set["literal"])

    return bnode_triples, ignored_triples


def get_fnode(s, o):
    """(8)Get the filled node of a subject and predicate pair, if one
    of which is a BNode.

    Parameters
    ----------
    s : rdflib.term.BNode|URIRef
    o : rdflib.term.BNode|Literal|URIRef

    Returns
    -------
    reflib.term.Literal|URIRef
    """
    if isinstance(s, BNode) and isinstance(o, BNode):
        raise Exception("Both s and o are blank")

    if not isinstance(s, BNode) and not isinstance(o, BNode):
        raise Exception("Both s and o are filled")

    if isinstance(s, BNode):
        return o

    else:
        return s


def load_triples_into_adb_graph(
    triples, adb_graph, vertex_collections, edge_collections, ro=None
):  # (9) <class 'arango.graph.Graph'>

    for s, p, o in triples:

        create_or_get_vertices_from_triple(
            adb_graph, vertex_collections, s, p, o, ro=ro
        )
        create_or_get_edge_from_triple(
            adb_graph, vertex_collections, edge_collections, s, p, o, ro=ro
        )

    for s, p, o in triples:

        update_vertex_from_triple(adb_graph, vertex_collections, s, p, o, ro=ro)


def create_or_get_vertex(
    adb_graph, vertex_collections, vertex_name, vertex_key, vertex_term
):
    """(10)Create, or get the identified vertex, creating vertex
    collections as needed.

    Parameters
    ----------
    adb_graph : arango.graph.Graph
        An ArangoDB graph instance
    vertex_collections : dict
        A dictionary with vertex name keys containing ArangoDB vertex
        collection instance values
    vertex_name : str
        The vertex collection name
    vertex_key : str
        The vertex key
    vertex_term : str
        The vertex ontology term

    Returns
    -------
    vertex : dict
        The ArangoDB vertex document
    """
    if vertex_name not in VALID_VERTICES:
        print(f"Skipping invalid vertex name: {vertex_name}")
        return

    vertex = {}

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

    else:
        vertex = vertex_collections[vertex_name].get(vertex_key)

    return vertex


def create_or_get_vertices_from_triple(
    adb_graph, vertex_collections, s, p, o, ro=None
):
    """(11)Create, or get vertices for the subject and object of the
    triple, creating vertex collections as needed.

    Parameters
    ----------
    adb_graph : arango.graph.Graph
        An ArangoDB graph instance
    vertex_collections : dict
        A dictionary with vertex name keys containing ArangoDB vertex
        collection instance values
    s : rdflib.term.BNode|URIRef
        Subject of triple
    p : rdflib.term.URIRef
        Predicate of triple
    o : rdflib.term.BNode|Literal|URIRef
        Object of triple

    Returns
    -------
    vertices : list(dict)
        List of ArangoDB vertex documents
    """
    if isinstance(o, Literal):
        # print(f"Skipping literal object in triple: {(s, p, o)}")
        return

    vertices = []

    for term in [s, o]:

        oid, number, term, _fragment, term_type = parse_term(term, ro=ro)

        if term_type != "class":
            continue

        vertex_name = oid
        vertex_key = number
        vertex_term = term

        vertex = create_or_get_vertex(
            adb_graph, vertex_collections, vertex_name, vertex_key, vertex_term
        )

        if vertex is None:
            # Message printed in previous function call
            return

        vertices.append(vertex)

    return vertices


def create_or_get_edge_from_triple(
    adb_graph, vertex_collections, edge_collections, s, p, o, ro=None
):  # (12)

    if isinstance(o, Literal):
        # print(f"Skipping literal object in triple: {(s, p, o)}")
        return

    s_oid, s_number, s_term, _s_fragment, s_term_type = parse_term(s, ro=ro)

    if s_term_type != "class":
        print(f"Skipping invalid subject type in triple: {(s, p, o)}")
        return

    from_vertex_name = s_oid
    from_vertex_key = s_number
    from_vertex_term = s_term

    _p_oid, _p_number, _p_term, p_fragment, p_term_type = parse_term(p, ro=ro)

    if not (
        p_term_type == "predicate"
        or (p_term_type == "class" and p_fragment is not None)
    ):
        print(f"Skipping invalid predicate type in triple: {(s, p, o)}")
        return

    predicate = p_fragment

    o_oid, o_number, o_term, _o_fragment, o_term_type = parse_term(o, ro=ro)

    if o_term_type != "class" and o_term_type != "literal":
        print(f"Skipping invalid object type in triple: {(s, p, o)}")
        return

    to_vertex_name = o_oid
    to_vertex_key = o_number
    to_vertex_term = o_term

    edge = create_or_get_edge(
        adb_graph,
        vertex_collections,
        edge_collections,
        from_vertex_name,
        from_vertex_key,
        from_vertex_term,
        to_vertex_name,
        to_vertex_key,
        to_vertex_term,
        predicate,
    )

    return edge


def create_or_get_edge(
    adb_graph,
    vertex_collections,
    edge_collections,
    from_vertex_name,
    from_vertex_key,
    from_vertex_term,
    to_vertex_name,
    to_vertex_key,
    to_vertex_term,
    predicate,
):  # (13) <class 'arango.collection.EdgeCollection'>

    from_vertex = create_or_get_vertex(
        adb_graph,
        vertex_collections,
        from_vertex_name,
        from_vertex_key,
        from_vertex_term,
    )

    if from_vertex is None:
        # Message printed in previous function call
        return

    to_vertex = create_or_get_vertex(
        adb_graph, vertex_collections, to_vertex_name, to_vertex_key, to_vertex_term
    )

    if to_vertex is None:
        # Message printed in previous function call
        return

    edge = {}

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

    else:
        edge = edge_collections[edge_name].get(edge_key)

    return edge


def update_vertex_from_triple(adb_graph, vertex_collections, s, p, o, ro=None):  # (14)

    if not isinstance(o, Literal):
        # print(f"Skipping non-literal object in triple: {(s, p, o)}")
        return

    s_oid, s_number, s_term, s_fragment, s_term_type = parse_term(s, ro=ro)

    if s_term_type != "class":
        print(f"Skipping invalid subject type in triple: {(s, p, o)}")
        return

    vertex_name = s_oid
    vertex_key = s_number
    vertex_term = s_term

    vertex = create_or_get_vertex(
        adb_graph, vertex_collections, vertex_name, vertex_key, vertex_term
    )

    if vertex is None:
        # Message printed in previous function call
        return

    p_oid, p_number, p_term, p_fragment, p_term_type = parse_term(p, ro=ro)

    if not (
        p_term_type == "predicate"
        or (p_term_type == "class" and p_fragment is not None)
    ):
        print(f"Skipping invalid predicate type in triple: {(s, p, o)}")
        return

    predicate = p_fragment

    if isinstance(o.value, datetime):
        value = str(o.value)

    else:
        value = o.value

    if not predicate in vertex:
        vertex[predicate] = value

    else:
        if not isinstance(vertex[predicate], list):
            vertex[predicate] = [vertex[predicate]]
        if value not in vertex[predicate]:
            vertex[predicate].append(value)

    vertex_collections[vertex_name].update(vertex)

    return vertex


def main():

    ols_dirname = (
        "/Users/raymondleclair/Projects/NLM/NLM-KN/springbok-nlm-kn/nlm-kn/data/ols"
    )

    bioportal_dirname = "/Users/raymondleclair/Projects/NLM/NLM-KN/springbok-nlm-kn/nlm-kn/data/bioportal"

    ontologies = [
        {
            "cl_filename": "general_cell_types_upper_slim.owl",
            "db_name": "BioPortal-BNode",
            "graph_name": "CL-BNode",
        },
        {
            "cl_filename": "general_cell_types_upper_slim.owl",
            "db_name": "BioPortal-Slim",
            "graph_name": "CL-Slim",
        },
        # {"cl_filename": "cl.owl", "db_name": "BioPortal-Full", "graph_name": "CL-Full"},
    ]

    for ontology in ontologies:

        cl_filename = ontology["cl_filename"]
        db_name = ontology["db_name"]
        graph_name = ontology["graph_name"]

        log_filename = f"{graph_name}.log"

        ro_filename = "ro.owl"
        ro, _ = parse_ols(ols_dirname, ro_filename)

        _, ids = parse_ols(bioportal_dirname, cl_filename)
        print(ids)

        rdf_graph = Graph()
        rdf_graph.parse(Path(bioportal_dirname) / cl_filename)

        triples = []
        triples_filename = log_filename.replace(".log", "_triples.txt")
        with open(triples_filename, "w") as fp:
            for triple in rdf_graph:
                triples.append(triple)
                fp.write(str(triple) + "\n")

        triple_types = count_triple_types(rdf_graph)
        pprint(triple_types)

        fnode_triples = collect_fnode_triples(rdf_graph)

        fnode_triples_filename = log_filename.replace(".log", "_fnode_triples.txt")
        with open(fnode_triples_filename, "w") as fp:
            for fnode_triple in fnode_triples:
                fp.write(str(fnode_triple) + "\n")

        bnode_triple_sets = {}

        collect_bnode_triple_sets(rdf_graph, bnode_triple_sets, use="subject", ro=ro)
        collect_bnode_triple_sets(rdf_graph, bnode_triple_sets, use="object", ro=ro)

        bnode_triple_sets_filename = log_filename.replace(
            ".log", "_bnode_triple_sets.txt"
        )
        with open(bnode_triple_sets_filename, "w") as fp:
            pprint(bnode_triple_sets, fp)

        bnode_triples, ignored_triples = create_bnode_triples_from_bnode_triple_sets(
            bnode_triple_sets, ro=ro
        )

        bnode_triples_filename = log_filename.replace(".log", "_bnode_triples.txt")
        with open(bnode_triples_filename, "w") as fp:
            for bnode_triple in bnode_triples:
                fp.write(str(bnode_triple) + "\n")

        fnode_triples.extend(bnode_triples)

        adb.delete_database(db_name)
        db = adb.create_or_get_database(db_name)

        adb.delete_graph(db, graph_name)
        adb_graph = adb.create_or_get_graph(db, graph_name)

        vertex_collections = {}
        edge_collections = {}

        if db_name == "BioPortal-BNode":
            VALID_VERTICES.extend(["BNode", "RO"])
            triples_to_populate = triples
        else:
            triples_to_populate = fnode_triples
        load_triples_into_adb_graph(
            triples_to_populate, adb_graph, vertex_collections, edge_collections, ro=ro
        )


if __name__ == "__main__":
    main()
