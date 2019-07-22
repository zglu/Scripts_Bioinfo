## Functional enrichment in given genes

These are scripts to perform functional enrichment in input genes

#### For Pfam, InterPro, SMART etc use:

    ./_enrichmentTest.sh [list of test genes] [list of background genes] [gene annotation] [annotation names]

all are sorted and one gene per line (the annotations should be combined to the same gene line)


#### For KEGG pathways:

First need to get a pathway list that suitable for your testing. The following command will get a list of all pathways in the current KEGG database, edit the list to fit your need.

    curl -# http://rest.kegg.jp/list/pathway > ./pathway.list

Next we need to get the K orthologues in each map

    mkdir maps
    ./getMaps.sh 

Finally perform enrichment test in genes of interest


    ./kegg_enrichmentTest.sh [list of test genes] [list of background genes] [gene-K annotation] [pathway names]


#### For GO terms use:

    Rscript _TopGO.R [list of test genes] [gene-GO annotation file]

Can change the node size, top genes etc in the script.
