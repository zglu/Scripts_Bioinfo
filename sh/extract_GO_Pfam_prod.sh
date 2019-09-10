#!/bin/bash

## extract GO terms, Pfam, and InterPro accessions and product information from Sm_v7.x_all_validated.gff
## ./extract_go_Pfam.sh [Sm_v7.x_all_validated.gff]
# GO
cat $1 | grep Ontology_term | awk -F "\t" '{print $9}' | sed 's/ /_/g' |sed 's/;/ /g'| awk '{print $2 " " $NF}'| sed 's/Parent=//; s/Ontology_term=//'| sort | sed 's/,/ /g'| awk -v OFS='\t' '{for (i=2;i<=NF;i++) print $1,$i}' | sort -u | awk '$1!=p{if(p)print s; p=$1; s=$0; next}{sub(p,x); s=s $0} END{print s}' > all_GO.txt

# InterPro
grep =InterPro $1 | awk -F "\t" '{print $9}' | sed 's/ /_/g' |sed 's/;/ /g'| awk '{print $2 " " $NF}'|grep InterPro |sed 's/Parent=//; s/polypeptide_domain=InterPro://' | sed 's/_Pfam/ /'|awk '{print $1 " " $2}'| sed 's/,/ /g'| awk -v OFS='\t' '{for (i=2;i<=NF;i++) print $1,$i}' > ip1.txt
grep =InterPro $1 | awk -F "\t" '{print $9}' | sed 's/ /_/g' |sed 's/;/ /g'| awk '{print $2 " " $(NF-1)}'|grep InterPro |sed 's/Parent=//; s/polypeptide_domain=InterPro://' | sed 's/_Pfam/ /'|awk '{print $1 " " $2}'| sed 's/,/ /g'| awk -v OFS='\t' '{for (i=2;i<=NF;i++) print $1,$i}' | cat - ip1.txt | sort -u | awk '$1!=p{if(p)print s; p=$1; s=$0; next}{sub(p,x); s=s $0} END{print s}' | tr '\t' ',' | sed 's/,/ /' > all_interpro.txt 
rm -f ip1.txt

# Pfam
cat $1 |grep 'Pfam:'|awk -F "\t" '{print $9}' | sed 's/ /_/g' |sed 's/;/ /g'| awk '{print $2 " " $NF}'|grep Pfam |sed 's/Parent=//; s/Pfam://'| sed 's/_PF/ PF/'|sort| awk '{print $1 " " $3}'| sort | sed 's/,/ /g'| awk -v OFS='\t' '{for (i=2;i<=NF;i++) print $1,$i}' | sort -u | awk '$1!=p{if(p)print s; p=$1; s=$0; next}{sub(p,x); s=s $0} END{print s}' > pf1.txt
cat $1 |grep 'Pfam:'|awk -F "\t" '{print $9}' | sed 's/ /_/g' |sed 's/;/ /g'| awk '{print $2 " " $(NF-1)}'|grep Pfam |sed 's/Parent=//; s/Pfam://'| sed 's/_PF/ PF/'|sort| awk '{print $1 " " $3}'| sort | sed 's/,/ /g'| awk -v OFS='\t' '{for (i=2;i<=NF;i++) print $1,$i}' | sort -u > pf2.txt
cat pf1.txt pf2.txt | sort -u | awk '$1!=p{if(p)print s; p=$1; s=$0; next}{sub(p,x); s=s $0} END{print s}' | sed 's/ /,/g'| sed 's/,/ /' > all_Pfam.txt
rm -f pf1.txt pf2.txt

echo "Unique GO terms:"
cat all_GO.txt | sed 's/,/ /g'| awk -v OFS='\t' '{for (i=2;i<=NF;i++) print $1,$i}' | awk '{print $2}'| sort -u | wc -l

echo "Unique Pfam terms:"
cat all_Pfam.txt | sed 's/,/ /g'| awk -v OFS='\t' '{for (i=2;i<=NF;i++) print $1,$i}' | awk '{print $2}'| sort -u | wc -l

## extract product information from _all_validated.gff
cat $1 | grep mRNA | awk -F '\t' '{print $9}'| sed 's/;/%/g; s/ /_/g; s/%/ /g'|awk '{print $1 " " $3}'| sed 's/ID=//;s/product=//'|sort > all_prod.txt
