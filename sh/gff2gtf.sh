#!/usr/local/bin/bash

# usage: ./gff2gtf.sh [GFF]
grep -v UTR $1 > $1-noUTR.gff
gffread $1-noUTR.gff -T -o $1.gtf

rm -f $1-noUTR.gff
