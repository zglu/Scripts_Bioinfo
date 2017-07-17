#! /usr/bin/env python

# zl3@sanger.ac.uk // last update: 13/07/2017
# add notes to the end of 'transcript' line in a gff3 file
# such as ....ID=Smp_100010.1;Parent=Smp_100010 --> ...Parent=Smp_100010;note=psid:Smp_12345
# USAGE: pyton ~zl3/scripts/Pair_Replace.py <id notes> <gff> <output>

import sys
import re

IDPAIRS = sys.argv[1] # tab separated file like: Smp_100010     Smp_12345
GFF = sys.argv[2] # gff file
FILEOUT = sys.argv[3] # output file

d = {}
with open(IDPAIRS) as f1:
    for line1 in f1:
        line1 = line1.rstrip()
        (aug, smp) = line1.split("\t")
        d['Parent=' + aug] = smp

with open(GFF) as f2:
    with open(FILEOUT, 'w') as fo:
        for line2 in f2:
            line2 = line2.rstrip()
            pat = re.compile(r'\b(' + '|'.join(d.keys()) + r')$') #\b(key1|key2|key...)$
            s = pat.sub(lambda x: x.group() + ';note=psid:' + d[x.group()], line2)   # replace 'Parent=Smp_100010' with 'Parent=g1;note=psid:Smp_12345'
            fo.write(s + "\n")
    fo.close()
