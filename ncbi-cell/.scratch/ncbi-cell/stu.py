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

        n_triples = 0
        n_triples_with_one_bnode = 0
        n_triples_with_double_bnodes = 0

        triple_sets = {}

        for s, p, o in v["graph"]:
            n_triples += 1

            if isinstance(s, BNode) and isinstance(o, BNode):
                n_triples_with_double_bnodes += 1
                continue

            if isinstance(s, BNode):
                if s not in triple_sets:
                    n_triples_with_one_bnode += 1
                    triple_sets[s] = {}
                    triple_sets[s]["relation"] = []
                    triple_sets[s]["annotation"] = []
                    triple_sets[s]["literal"] = []
                    triple_sets[s]["class"] = []
                    triple_sets[s]["other"] = []

                _, _, _, _, s_term_type = parse_term(s)
                _, _, _, fragment, _ = parse_term(p)
                _, _, _, _, o_term_type = parse_term(o)

                if fragment in ["someValuesFrom", "onProperty", "subClassOf"]:
                    triple_sets[s]["relation"].append((s, p, o))

                elif fragment in ["annotatedProperty", "annotatedTarget", "annotatedSource"]:
                    triple_sets[s]["annotation"].append((s, p, o))

                elif fragment in ["hasDbXref", "source"]:
                    triple_sets[s]["literal"].append((s, p, o))

                elif s_term_type == "class" or o_term_type == "class":
                    triple_sets[s]["class"].append((s, p, o))

                else:
                    triple_sets[s]["other"].append((s, p, o))

        for s, p, o in v["graph"]:

            if isinstance(s, BNode) and isinstance(o, BNode):
                continue

            if isinstance(o, BNode):
                if o not in triple_sets:
                    n_triples_with_one_bnode += 1
                    triple_sets[o] = {}
                    triple_sets[o]["relation"] = []
                    triple_sets[o]["annotation"] = []
                    triple_sets[o]["literal"] = []
                    triple_sets[o]["class"] = []
                    triple_sets[o]["other"] = []

                _, _, _, _, s_term_type = parse_term(s)
                _, _, _, fragment, _ = parse_term(p)
                _, _, _, _, o_term_type = parse_term(o)

                if fragment in ["someValuesFrom", "onProperty", "subClassOf"]:
                    triple_sets[o]["relation"].append((s, p, o))

                elif fragment in ["annotatedProperty", "annotatedTarget", "annotatedSource"]:
                    triple_sets[o]["annotation"].append((s, p, o))

                elif fragment in ["hasDbXref", "source"]:
                    triple_sets[o]["literal"].append((s, p, o))

                elif s_term_type == "class" or o_term_type == "class":
                    triple_sets[o]["class"].append((s, p, o))

                else:
                    triple_sets[o]["other"].append((s, p, o))

        pprint(triple_sets)

        print(f"n_triples: {n_triples}")
        print(f"n_triples_with_one_bnode: {n_triples_with_one_bnode}")
        print(f"n_triples_with_double_bnodes: {n_triples_with_double_bnodes}")
