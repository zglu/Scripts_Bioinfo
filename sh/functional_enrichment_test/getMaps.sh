#!/bin/bash

# get the list of all pathways
# (curl -# http://rest.kegg.jp/list/pathway > ./pathway.list)
# remove path: in the list to have a valid name for the following command
# remove global and unrelated maps:
# 011 global map (lines linked to KOs)
# 012 overview map (lines linked to KOs)
# 010 chemical structure map (no KO expansion)
# 07 drug structure map (no KO expansion)
#
# 05 Human Diseases

if [ $? -eq 0 ]; then
    for next in $(cut -f1 pathway.list); do
    (curl -# http://rest.kegg.jp/link/ko/$next > maps/$next-K.txt)
    if [ $? -eq 0 ]; then
        echo "Retrieved $next ko map"   
    else
        echo "There was a problem in data retrieval"
        exit 1
    fi
    done
    exit 0
else 
    echo "There was a problem in data retrieval"
    exit 1
fi
