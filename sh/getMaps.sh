#!/bin/bash
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
