#! /usr/bin/env python

# zl3@sanger.ac.uk // last update: 13/07/2017
# replace the whole word in file2 (I use a gff3 file) with the paird word in file1 (id pairs file)
# e.g., replaces g10 to Smp_100010 and g100 to Smp_100100; and g10.t1 to g10.1
# USAGE: pyton ~zl3/scripts/Pair_Replace.py <id-pairs> <gff> <output>

import sys
import re

IDPAIRS = sys.argv[1] # tab separated file
GFF = sys.argv[2] # gff file
FILEOUT = sys.argv[3] # output file

d = {}
dt = {}
with open(IDPAIRS) as f1:
    for line1 in f1:
        line1 = line1.rstrip()
        (aug, smp) = line1.split("\t")
        d[aug] = smp
        dt[(aug + '.t')] = smp

with open(GFF) as f2:
    with open(FILEOUT, 'w') as fo:
        for line2 in f2:
            line2 = line2.rstrip()
            pat1 = re.compile(r'\b(' + '|'.join(dt.keys()) + r')') # \b(key1.t|key2.t|key...)
            pat2 = re.compile(r'\b(' + '|'.join(d.keys()) + r')\b') # \b(key1|key2|key...)\b
            s1 = pat1.sub(lambda x: dt[x.group()] + '.', line2) # replace 'g1.t1' with 'Smp_1.1'
            s2 = pat2.sub(lambda x: d[x.group()], s1)   # replace 'g1' with 'Smp_1'
            fo.write(s2 + "\n")
    fo.close()
