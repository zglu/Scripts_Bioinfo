#!/bin/bash

## split a fasta file with multiple sequencs to single sequences
## split all: ./splitfa.sh input.fa
for id in $(grep '>' $1 | cut -d \> -f2); do
## or partial: ./splitfa.sh input.fa ids
# for id in $(cut -f1 $2); do
    (cat $1 |sed /$id/',/>/!d'| awk '/>/{i++}i==1' > $id.fa)
done
