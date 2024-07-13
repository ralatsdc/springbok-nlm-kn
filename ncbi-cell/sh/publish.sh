#!/usr/bin/env bash

org_files=$(ls *.org)
for org_file in $org_files; do
    ipynb_file=$(echo $org_file | sed s/.org/.ipynb/)

    echo "Tangling $org_file"
    emacs --batch --eval "(require 'org)" --eval "(org-babel-tangle-file \"$org_file\")"

    echo "Converting $org_file to $ipynb_file"
    cat $org_file | sed 's/src python/src jupyter-python/' > tmp.org
    pandoc tmp.org -o $ipynb_file
    rm tmp.org
    
done

black *.py
