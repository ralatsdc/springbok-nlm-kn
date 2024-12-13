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
        mid = annot.loc[marker_name, "ensembl_gene_id"]
        if isinstance(mid, pd.core.series.Series):
            marker_ids[marker_name] = mid.to_list()
        else:
            marker_ids[marker_name] = [mid]
    else:
        print(f"Could not find marker name: {marker_name}")

# Get resources for each marker name
resources = [
    "diseases",
    "drugs",
    "interactions",
    "tractability",
    "expression",
    "depmap",
]
gdata = {}
gdata["resources"] = resources
# nMNm = 0
for mname, mids in marker_ids.items():
    # nMNm += 1
    gdata[mname] = {}
    gdata[mname]["ids"] = mids

    # Get resources for each marker id
    for mid in mids:
        gdata[mname][mid] = {}
        for resource in resources:
            try:
                gdata[mname][mid][resource] = gget.opentargets(
                    mid, resource=resource, json=True, verbose=True
                )
            except Exception as exc:
                gdata[mname][mid][resource] = {}
                print(
                    f"Could not gget resource: {resource} for marker id: {mid} for marker name: {mname}"
                )

    # if nMNm == 8:
    #     break

# Write the resources for all marker names to a JSON file
gdata_filename = mdata_filename.replace(".xlsm", ".json")
with open(mdata_dirname / gdata_filename, "w") as ofp:
    json.dump(gdata, ofp, indent=4)
