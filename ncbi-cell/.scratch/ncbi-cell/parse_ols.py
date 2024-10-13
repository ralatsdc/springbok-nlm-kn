from pathlib import Path

from lxml import etree

from jkl import parse_term

owl = "{http://www.w3.org/2002/07/owl#}"
rdf = "{http://www.w3.org/1999/02/22-rdf-syntax-ns#}"
rdfs = "{http://www.w3.org/2000/01/rdf-schema#}"

ols_dir = Path(
    "/Users/raymondleclair/Projects/NLM/NLM-KB/springbok-ncbi-cell/ncbi-cell/data/ols"
)

ols_fnm = "iao.owl"
ols_fnm = "ro.owl"
ols_fnm = "bfo_classes_only.owl"

root = etree.parse(str(ols_dir / ols_fnm))

ols = {}

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

            ols[term] = label
