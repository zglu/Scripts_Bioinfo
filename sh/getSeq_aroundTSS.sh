#!/bin/bash

## get the up- and down-stream sequences from TSS
## usage: ./getSeq_aroundTSS.sh [genes] [upstream region] [downstream region]

grep -f $1 ~/Documents/SCHISTO/V7/annotations/Smv7.3/Sm_v7.3_Basic-noGaps.gff > $1.gff
python3 ~/Google\ Drive/Scripts_Bioinfo/py/up-downstream_sequence.py $1.gff ~/Documents/SCHISTO/V7/Smansoni_v7.fa gene $2 $3 > ${1}_up$2-down$3.fa
