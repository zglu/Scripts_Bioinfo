#!/bin/bash

# get upstream and downstream regions around TSS
# usage: ./aroundTSS_Bed.sh [upstream_region] [downstream_region]
declare -i UPREGION=$1
declare -i DWREGION=$2

cat Smv7.2.genes.bed | awk '{if ($6=="+") print $1, $2-"'$UPREGION'", $2+"'$DWREGION'", $4, $5, $6; else print $1, $3-"'$DWREGION'", $3+"'$UPREGION'", $4, $5, $6}'| awk '{if($2<0) print $1, "1", $3, $4, $5, $6; else print $0}'|tr ' ' '\t' > promoter_up"$UPREGION"-down"$DWREGION".bed
