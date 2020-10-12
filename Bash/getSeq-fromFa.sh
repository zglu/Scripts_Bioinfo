#!/bin/bash
    for id in $(cut -f1 s7newids); do
    (cat WBPS7.Sm.protein.fa |sed /$id/',/>/!d'| awk '/>/{i++}i==1' > ./$id.fa)
    done
