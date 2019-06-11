#!/bin/bash

### Perform function enrichment (Pfam, InterPro) analysis via hypergeometric test
### Need to provide sorted files: test genes ($1), background genes ($2), gene annotation ($3), term names ($4) 

# gene annotation
# Smp_000020 PF07555
# Smp_000040 PF13374,PF13424

# term names (two columns only)
# PF00007	Cys_knot:Cystine-knot_domain
# PF00008	EGF:EGF-like_domain


## get feasible number of genes for test and background
TESTFEA="$(join $1 $3 | wc -l | awk '{print $1}')"
ALLFEA="$(join $2 $3 | wc -l | awk '{print $1}')"

## get list of function domains and number of mapped genes
join $2 $3 | sed 's/,/ /g'| awk -v OFS='\t' '{for (i=2;i<=NF;i++) print $1,$i}' | awk '{print $2}'| sort | uniq -c  | awk '{print $2 " " $1}'> $2-mapped

join $1 $3 | sed 's/,/ /g'| awk -v OFS='\t' '{for (i=2;i<=NF;i++) print $1,$i}' | awk '{print $2}'| sort | uniq -c  | awk '{print $2 " " $1}' > tmp1.txt
join $1 $3 | sed 's/,/ /g'| awk -v OFS='\t' '{for (i=2;i<=NF;i++) print $1,$i}' | awk '{print $2 "\t" $1}'| sort | awk '$1!=p{if(p)print s; p=$1; s=$0; next}{sub(p,x); s=s $0} END{print s}' | tr '\t' ','| sed 's/,/ /'| join tmp1.txt - > $1-mapped

join $1-mapped $2-mapped | awk '{print $1, $2, $4, "'$TESTFEA'"-$2, "'$ALLFEA'"-$4, $3}' | awk '$2>0' | join $4 - > $1-fisherTable.txt

rm -f tmp1.txt

# run hypergeometric test
Rscript _FisherTest.R $1

