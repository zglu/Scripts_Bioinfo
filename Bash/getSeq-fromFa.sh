#!/bin/bash

# usage: ./getSeq-fromFa.sh [id] [fasta file]

sed /$1/',/>/!d' $2 | awk '/>/{i++}i==1'

#eg: sed /Smp_076300/',/>/!d' ~/zl3/SCHISTO/v9/Sm_v9.pep.fa | awk '/>/{i++}i==1'

# batch
#for id in $(cut -f1 IDLIST); do
#(sed /$id/',/>/!d' FASTAFILE| awk '/>/{i++}i==1' > ./$id.fa)
#done

