#!/bin/zsh

# Check for results from tutorial
if [ ! -f results-from-tutorial/cluster_results.csv ]; then
    echo "Results from tutorial not found"
    echo "Exiting ..."
    exit 1
else
    echo "Results from tutorial found"
fi

# Check for results from test
if [ ! -f results-from-test/cluster_results.csv ]; then
    echo "Results from test not found"
    echo "Running NS-Forest ..."
    source ../../.venv/bin/activate
    ./nsforest.py --run-nsforest-with-preprocessing -c cluster -d results-from-test/ ./data-for-test/adata_layer1.h5ad
else 
   echo "Results from test found"
fi

# Compare results
if diff results-from-tutorial/cluster_results.csv results-from-test/cluster_results.csv ; then
   echo "PASS: Tutorial and test results identical"
else
   echo "FAIL: Tutorial and test results different"
fi
