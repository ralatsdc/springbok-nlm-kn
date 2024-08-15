#!/usr/bin/env python

import json
from pathlib import Path
import re


def write_markdown(fp, cell):
    """Write a Markdown cell.

    Parameters
    ----------
    fp : _io.TextIOWrapper
        File pointer
    cell : dict
        The Markdown cell dictionary

    Returns
    -------
    None
    """
    # Write each line
    spcs = ""
    in_lisp = False
    in_bash = False
    for line in cell["source"]:

        # Substitute paired double asterisks with paired single
        # asterisks
        m = re.search(r"(\*\*).*(\*\*)", line)
        if m:
            line = re.sub(r"(\*\*).*(\*\*)", m.group(0).replace("**", "*"), line)

        # Remove backslash escapes from paired left and right angle
        # brackets
        m = re.search(r"(\\\<).*(\\\>)", line)
        if m:
            line = re.sub(r"(\\\<).*(\\\>)", m.group(0).replace("\\", ""), line)

        # Replace leading sharps with asterisks
        m = re.match(r"(#+)(?!\+)", line)
        if m:
            line = re.sub(r"(^#+)(?!\+)", m.group(1).replace("#", "*"), line)

        # Replace Markdown with Emacs Org mode link format. Note that
        # order of link and description is reversed, and "file"
        # protocol needs to be added.
        m = re.search(r"(\[.*\])(\(.*?\))", line)
        if m:
            link = m.group(2).replace("(", "[").replace(")", "]")
            if not re.match(r"\[http", link):
                link = link.replace("[", "[file:")
            desc = m.group(1)
            line = re.sub(
                r"(\[.*\])(\(.*?\))",
                f"[{link}{desc}]",
                line,
            )

        # Replace three spaces after leading hyphen with one
        line = re.sub(r"^-   ", "- ", line)

        # Replace four leading spaces with two
        line = re.sub(r"^    ", "  ", line)

        # Handle emacs-lisp source blocks
        if "``` commonlisp" in line:
            in_lisp = True
            fp.write("#+begin_src emacs-lisp :session shared :results silent\n")
            spcs = "  "
            continue
        if in_lisp and "```" in line:
            in_lisp = False
            fp.write("#+end_src\n")
            spcs = ""
            continue

        # Handle shell source blocks
        if "``` bash" in line:
            in_bash = True
            fp.write("#+begin_src sh\n")
            spcs = "  "
            continue
        if in_bash and "```" in line:
            in_bash = False
            fp.write("#+end_src\n")
            spcs = ""
            continue

        # Substitute backticks with tildes
        line = line.replace("`", "~")

        fp.write(spcs + line)

    fp.write("\n\n")


def write_code(fp, cell):
    """Write a code cell.

    Parameters
    ----------
    fp : _io.TextIOWrapper
        File pointer
    cell : dict
        The code cell dictionary

    Returns
    -------
    None
    """
    # Create, then write, begin source line using metadata to populate
    # options
    begin_src_str = "#+begin_src python"
    for k, v in cell["metadata"].items():
        begin_src_str += f" :{k} {v}"
    fp.write(begin_src_str + "\n")

    # Write source lines
    spcs = "  "
    for line in cell["source"]:
        if re.fullmatch(r"\s*", line):
            fp.write(line)
        else:
            fp.write(spcs + line)

    # Write end source line
    fp.write("#+end_src\n\n")


def publish_ipynb_to_org(ipynb_path):
    """Load the specified Jupyter notebook, and write the contents as
    an Org mode file.

    Parameters
    ----------
    ipynb_filename : Path
       Path of the Jupyter notebook file

    Returns
    -------
    None
    """
    # Load Jupyter notebook file
    with open(ipynb_path, "r") as fp:
        data = json.load(fp)

    # Write Org mode file
    org_path = ipynb_path.parents[1] / "org" / ipynb_path.name.replace(".ipynb", ".org")
    print(f"Converting {ipynb_path.name} to {org_path.name}")
    with open(org_path, "w") as fp:
        for cell in data["cells"]:
            if cell["cell_type"] == "markdown":
                write_markdown(fp, cell)
            elif cell["cell_type"] == "code":
                write_code(fp, cell)


def main():
    """Publish each Jupyter notebook in the directory containing this
    script to its corresponding Org mode file.

    Parameters
    ----------
    None

    Returns
    -------
    None
    """
    ipynb_home = Path(__file__).parent.resolve()
    for ipynb_path in ipynb_home.glob("*.ipynb"):
        publish_ipynb_to_org(ipynb_path)


if __name__ == "__main__":
    main()
