#!/usr/local/bin/bash
    for id in $(cut -d \| -f2 Ratt_AugSTAR/SmpAug_ids.txt); do
    (cat ../Aug_STAR_all.pep.fasta | sed /$id/',/>/!d'| awk '/>/{i++}i==1' > Ratt_AugSTAR/Aug-fa/$id.fa)
    done
