#!/usr/bin/env bash

# Assign home directory of this script
home=$(dirname $(dirname $(realpath $0)))

# Tangle, and convert to .ipynb eac Org mode file
org_files=$(ls $home/org/*.org)
for org_file in $org_files; do
    ipynb_file=$(basename $org_file | sed s/.org/.ipynb/)

    echo "Tangling $org_file"
    emacs --batch --eval "(require 'org)" --eval "(org-babel-tangle-file \"$org_file\")"

    echo "Converting $org_file to $ipynb_file"
    cat $org_file | sed 's/src python/src jupyter-python/' > tmp.org
    pandoc tmp.org -o $home/ipynb/$ipynb_file
    rm tmp.org

done

# Blacken all Python files
black $home/py/*.py
