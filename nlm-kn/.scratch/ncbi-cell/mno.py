from pathlib import Path
from pprint import pprint

from rdflib import Graph
from rdflib.term import BNode, Literal, URIRef

from jkl import parse_term, PREDICATE_CLASSES

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

        triple_sets = {}
        bnodes = []
        for s, p, o in v["graph"]:

            if isinstance(s, BNode):
                if s not in triple_sets:
                    triple_sets[s] = []
                triple_sets[s].append((s, p, o))
                if o == URIRef("http://purl.obolibrary.org/obo/CL_0000235"):
                    bnodes.append(s)

            if isinstance(o, BNode):
                if o not in triple_sets:
                    triple_sets[o] = []
                triple_sets[o].append((s, p, o))
                if s == URIRef("http://purl.obolibrary.org/obo/CL_0000235"):
                    bnodes.append(o)

        class_triple_sets = {}
        for bnode, triple_set in triple_sets.items():
            for s, p, o in triple_set:
                if isinstance(s, BNode):
                    term_tuple = parse_term(o)

                else:
                    term_tuple = parse_term(s)

                if len(term_tuple) == 3:
                    if bnode not in class_triple_sets:
                        class_triple_sets[bnode] = []
                    class_triple_sets[bnode].append((s, p, o))

        class_relation_triple_sets = {}
        for bnode, class_triple_set in class_triple_sets.items():
            if len(class_triple_set) == 3:
                for s, p, o in class_triple_set:
                    if isinstance(s, BNode):
                        term_tuple = parse_term(o)

                    else:
                        term_tuple = parse_term(s)

                    if term_tuple[0] == "RO":
                        class_relation_triple_sets[bnode] = class_triple_set
                        break

        pprint(class_relation_triple_sets)

        for bnode, class_relation_triple_set in class_relation_triple_sets.items():

            s_prime = BNode()
            p_prime = BNode()
            o_prime = BNode()

            for s, p, o in class_relation_triple_set:
                if isinstance(s, BNode):
                    c = o

                else:
                    c = s

                term_tuple = parse_term(c)

                if term_tuple[0] == "CL":
                    s_prime = c

                if term_tuple[0] == "RO":
                    p_prime = c

                else:
                    o_prime = c

            print((s_prime, p_prime, o_prime))
