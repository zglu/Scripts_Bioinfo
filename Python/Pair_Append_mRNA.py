#! /usr/bin/env python

# zl3@sanger.ac.uk // last update: 13/07/2017
# add psid notes to the end of 'mRNA' line in a gff3 file
# such as ....ID=Smp_100010.1;Parent=Smp_100010 --> ID=Smp_100010;Parent=Smp_100010;product=Hypothetical protein...
# USAGE: python ~zl3/scripts/Pair_Replace.py <id notes> <gff> <output>

import sys
import re

IDPAIRS = sys.argv[1] # tab separated file like: Smp_300010.1	hypothetical protein
GFF = sys.argv[2] # gff file
FILEOUT = sys.argv[3] # output file

d = {}
with open(IDPAIRS) as f1:
    for line1 in f1:
        line1 = line1.rstrip()
        (smp, note) = line1.split("\t")
        d['ID=' + smp + ';'] = note # mRNA ID as key: eg. ID=Smp_100010.1;

with open(GFF) as f2:
    with open(FILEOUT, 'w') as fo:
        for line2 in f2:
            line2 = line2.rstrip()
            pat = re.compile(r'\b(' + '|'.join(d.keys()) + r')' + 'Parent=' + r'\d+$') #\b(key1|key2|key...)$
            s = pat.sub(lambda x: x.group() + ';product=' + d[x.group()], line2)   # replace ID=Smp_100010.1;Parent=Smp_100010 --> ID=Smp_100010;Parent=Smp_100010;product=Hypothetical protein...
            fo.write(s + "\n")
    fo.close()
