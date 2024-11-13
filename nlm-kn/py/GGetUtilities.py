import ast
import json
from pathlib import Path

import gget
import pandas as pd
import scanpy as sc

# Retrieve gene symbols
annot = sc.queries.biomart_annotations(
    "hsapiens", ["ensembl_gene_id", "external_gene_name"], use_cache=True
).set_index("external_gene_name")

# Load the manually curated data
mdata_dirname = Path("../data/nlm-kn")
mdata_filename = "HLCA_CellRef_matching_ver3_import1.xlsm"
mdata = pd.read_excel(mdata_dirname / mdata_filename, header=1, skiprows=0)

# Create a list of marker names
marker_names = []
[
    marker_names.extend(mns)
    for mns in mdata["HLCA_NSForestMarkers"].dropna().apply(ast.literal_eval)
]

# Create a dictionary of marker ids, using marker names as keys
marker_ids = {}
for marker_name in marker_names:
    if marker_name in annot.index:
        marker_id = annot.loc[marker_name, "ensembl_gene_id"]
        if isinstance(marker_id, pd.core.series.Series):
            # TODO: Decide how to handle multiple ids per name
            marker_id = marker_id.iloc[0]
        marker_ids[marker_name] = marker_id
    else:
        print(f"Could not find marker name: {marker_name}")

# Get diseases and drugs for each marker id
gdata = {}
for marker_name, marker_id in marker_ids.items():
    gdata[marker_name] = {}
    gdata[marker_name]["marker_id"] = marker_id
    try:
        gdata[marker_name]["diseases"] = gget.opentargets(
            marker_id, resource="diseases", json=True, verbose=True
        )
        gdata[marker_name]["drugs"] = gget.opentargets(
            marker_id, resource="drugs", json=True, verbose=True
        )
    except Exception as exc:
        print(f"Could not gget marker id: {marker_id}")

# Write the diseases and drugs data to a JSON file
gdata_filename = mdata_filename.replace(".xlsm", ".json")
with open(mdata_dirname / gdata_filename, "w") as ofp:
    json.dump(gdata, ofp)
