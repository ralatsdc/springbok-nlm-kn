from pathlib import Path
from urllib.parse import urlparse

from lxml import etree
import requests

owl = "{http://www.w3.org/2002/07/owl#}"

obo_dirpath = Path("../data/obo")
obo_purls = [
    "http://purl.obolibrary.org/obo/cl.owl",
    "http://purl.obolibrary.org/obo/ro.owl",
]

for obo_purl in obo_purls:

    obo_stem = Path(urlparse(obo_purl).path).stem
    obo_suffix = Path(urlparse(obo_purl).path).suffix

    print(f"Getting {obo_purl}")
    r = requests.get(obo_purl)

    obo_filepath_new = obo_dirpath / (obo_stem + "-new" + obo_suffix)

    print(f"Writing {obo_filepath_new}")
    with open(obo_filepath_new, "wb") as f:
        f.write(r.content)

    print(f"Parsing {obo_filepath_new}")
    root = etree.parse(obo_filepath_new)
    version_new = root.find(f"{owl}Ontology/{owl}versionInfo").text
    print(f"Found new version {version_new}")

    obo_filepath_cur = obo_dirpath / (obo_stem + obo_suffix)
    if obo_filepath_cur.exists():
        print(f"Parsing {obo_filepath_cur}")
        root = etree.parse(obo_filepath_cur)
        version_cur = root.find(f"{owl}Ontology/{owl}versionInfo").text
        print(f"Found current version {version_cur}")

        if version_new > version_cur:
            obo_filepath_old = obo_dirpath / (obo_stem + "-" + version_cur + obo_suffix)

            print(f"Renaming {obo_filepath_cur} to {obo_filepath_old}")
            obo_filepath_cur.rename(obo_filepath_old)

            print(f"Renaming {obo_filepath_new} to {obo_filepath_cur}")
            obo_filepath_new.rename(obo_filepath_cur)

        else:
            print(f"New version is not newer than current version")
            print(f"Removing {obo_filepath_new}")
            obo_filepath_new.unlink()

    else:
        print(f"Renaming {obo_filepath_new} to {obo_filepath_cur}")
        obo_filepath_new.rename(obo_filepath_cur)
