#!/bin/bash

## usage: ./kegg_enrichmentTest.sh [test genes] [background genes] [gene-K annotation] [pathway names]

### Perform function enrichment (KEGG pathways) analysis via hypergeometric test
### Need to provide sorted files: test genes ($1), background genes ($2), gene-K annotation ($3), pathway names ($4)

# gene-K annotation (kegg-7.2.txt)
# Smp_000020 K15719
# Smp_000040 K10407

# pathway names (two columns only; keggRef.txt)
# map00010	Glycolysis_/_Gluconeogenesis
# map00020	Citrate_cycle_(TCA_cycle)

### Need to have a folder maps/ with all K orthologs in each map
## - pathway.list contains pathways that you would like to find enrichment in
## getMaps.sh for retrieving K in each map

# map00010-K.txt
# path:map00010	ko:K00001
# path:map00010	ko:K00002
# path:map00010	ko:K00016

join $1 $3 | grep K | awk '{print $2 " " $1}'| sort > $1-K
join $2 $3 | grep K | awk '{print $2 " " $1}'| sort > $2-K

TESTMAPPED="$(wc -l $1-K | awk '{print $1}')"
BKGDMAPPED="$(wc -l $2-K | awk '{print $1}')"

## get gene counts in each map for test and background
for i in maps/*-K.txt; do printf $i; cat $i | sed 's/:/ /g' | awk '{print $4}' | join - $2-K | wc -l; done | sed 's/maps\///' > $2.mapped.counts
for i in maps/*-K.txt; do printf $i; cat $i | sed 's/:/ /g' | awk '{print $4}' | join - $1-K | wc -l; done | sed 's/maps\///' > $1.mapped.counts.par1
## get genes for each map in the test set
for i in maps/*-K.txt; do printf $i "\t"; cat $i | sed 's/:/ /g' | awk '{print $4}' | join - $1-K | awk '{print "map" " " $2}' | awk '$1!=p{if(p)print s; p=$1; s=$0; next}{sub(p,x); s=s $0} END{print s}'| tr ' ' ',' | sed 's/,/ /'| awk '{print " " $2}'; done | sed 's/maps\///' > $1.mapped.counts.par2
join $1.mapped.counts.par1 $1.mapped.counts.par2 > $1.mapped.counts
rm -f $1.mapped.counts.par1 $1.mapped.counts.par2

## construct the test table
join $1.mapped.counts $2.mapped.counts | sed 's/-K\.txt//'| awk '{print $1 " " $2 " " $4 " " "'$TESTMAPPED'"-$2 " " "'$BKGDMAPPED'"-$4 " " $3}' | awk '$2>0'|sort | join $4 - > $1-fisherTable.txt

# run hypergeometric test
Rscript _FisherTest.R $1

rm -f $1-K $2-K $1.mapped.counts $2.mapped.counts $1-fisherTable.txt
