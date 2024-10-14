from pathlib import Path

from lxml import etree

from jkl import parse_term

base = "{http://purl.obolibrary.org/obo/iao.owl}"
dc = "{http://purl.org/dc/elements/1.1/}"
obo = "{http://purl.obolibrary.org/obo/}"
owl = "{http://www.w3.org/2002/07/owl#}"
rdf = "{http://www.w3.org/1999/02/22-rdf-syntax-ns#}"
xml = "{http://www.w3.org/XML/1998/namespace}"
xsd = "{http://www.w3.org/2001/XMLSchema#}"
foaf = "{http://xmlns.com/foaf/0.1/}"
rdfs = "{http://www.w3.org/2000/01/rdf-schema#}"
swrl = "{http://www.w3.org/2003/11/swrl#}"
swrlb = "{http://www.w3.org/2003/11/swrlb#}"
terms = "{http://purl.org/dc/terms/}"
protege = "{http://protege.stanford.edu/plugins/owl/protege#}"
oboInOwl = "{http://www.geneontology.org/formats/oboInOwl#}"

ols_dir = Path(
    "/Users/raymondleclair/Projects/NLM/NLM-KB/springbok-ncbi-cell/ncbi-cell/data/ols"
)

iao_fnm = "iao.owl"

root = etree.parse(str(ols_dir / iao_fnm))

iao = {}

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
        if label_element is not None:
            label = label_element.text

            iao[term] = label