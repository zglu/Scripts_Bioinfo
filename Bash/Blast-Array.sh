#!/usr/local/bin/bash
declare -A idpairs # associative array

while IFS=\| read smp aug  # columns are seperated by pipe |
do
    idpairs[$aug]=$smp
done < SmpAug_ids.txt  # this is the id pair file

for aug in ${!idpairs[*]}
do  # make directories with single fa files, and another for the blast output
    blastp -query ../RattPacbio-fa/${idpairs[$aug]}.fa -subject Aug-fa/$aug.fa -out blast-out/${idpairs[$aug]}-$aug.txt
done
