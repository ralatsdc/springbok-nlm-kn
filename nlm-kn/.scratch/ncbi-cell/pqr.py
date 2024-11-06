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

        n_nodes_with_multiple_bnodes = 0
        for s, p, o in v["graph"]:

            if isinstance(s, BNode):

                if s not in triple_sets:
                    triple_sets[s] = {}
                    triple_sets[s]["BNodes"] = set()
                    triple_sets[s]["triples"] = []

                triple_sets[s]["BNodes"].add(s)
                if isinstance(o, BNode):
                    triple_sets[s]["BNodes"].add(o)

                if len(triple_sets[s]["BNodes"]) > 1:
                    n_nodes_with_multiple_bnodes += 1

                triple_sets[s]["triples"].append((s, p, o))


        n_object_nodes_not_in_triple_sets = 0
        for s, p, o in v["graph"]:

            if isinstance(o, BNode):

                if o not in triple_sets:
                    n_object_nodes_not_in_triple_sets += 1

    print(f"n_nodes_with_multiple_bnodes: {n_nodes_with_multiple_bnodes}")
    print(f"n_object_nodes_not_in_triple_sets: {n_object_nodes_not_in_triple_sets}")

    n_found = 0

    for bnode, triple_set in triple_sets.items():
        if len(triple_set["triples"]) < 3:
            continue
        
        n_classes_found = 0

        for s, p, o in triple_set["triples"]:
            if isinstance(s, BNode) and isinstance(o, BNode):
                n = p

            elif not isinstance(s, BNode):
                n = s
            
            else:  # Object must not be a Bnode
                n = o

            oid, number, term, fragment, term_type = parse_term(n)


            if term_type == "class":
                n_classes_found += 1
                
        if n_classes_found == 2:

            n_found += 1

            print("-----------------------------------")
            print(bnode, len(triple_set["triples"]))
            print("-----------------------------------")

            for s, p, o in triple_set["triples"]:
                print(parse_term(s), s)
                print(parse_term(p), p)
                print(parse_term(o), o)

    print(f"n_found: {n_found}")
