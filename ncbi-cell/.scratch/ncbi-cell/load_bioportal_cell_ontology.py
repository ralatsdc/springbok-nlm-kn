from pathlib import Path

from rdflib import Graph
from arango_rdf import ArangoRDF

import ArangoDB as adb


bioportal_dir = Path("/Users/raymondleclair/Projects/NLM/NLM-KB/springbok-ncbi-cell/ncbi-cell/data/bioportal")
bioportal_fnm = "general_cell_types_upper_slim.owl"
# bioportal_fnm = "cl.owl"
# bioportal_fnm = "owlapi.xrdf"
# bioportal_fnm = "pizza.owl"

rdf_graph = Graph()
rdf_graph.parse(str(bioportal_dir / bioportal_fnm))

adb.delete_database("BioPortal")
db = adb.create_or_get_database("BioPortal")

adbrdf = ArangoRDF(db)

adbrdf.rdf_to_arangodb_by_pgt(name="Cell Ontology", rdf_graph=rdf_graph, overwrite_graph=True)
