#! /usr/bin/env python

'''
Calculate UTR from gff containing ONLY mRNA and CDS
for Apollo gff files: coordinates are ordered
'''

# last update: 7/11/2017

import re
import sys

GFF = sys.argv[1] # single gene gff

with open(GFF) as fin:
    for line in fin:
        line = line.rstrip()
        line = line.split()
        lastCDS = None
        if "mRNA" in line:
            mStart = line[3]
            mEnd = line[4]
            i = 0
        elif "CDS" in line:
            lastCDS = line
            i += 1
            if i == 1:
                pStart = line[3]
                pStart_rv = line[4]
    line[8] = re.sub(";ID=(.+?)$", "", line[8])
# here lastCDS is the last line containing CDS
# if mRNA start/end = start_codon/stop_codon then may end up with -1
    if line[6] == "+":
        five_utr = [line[0], line[1], "five_prime_UTR", str(mStart), str(int(pStart) - 1), line[5], line[6], ".", line[8]]
        three_utr = [line[0], line[1], "three_prime_UTR", str(int(lastCDS[4]) + 1), str(mEnd), line[5], line[6], ".", line[8]]
    elif line[6] == "-":
        five_utr = [line[0], line[1], "five_prime_UTR", str(int(pStart_rv) + 1), str(mEnd), line[5], line[6], ".", line[8]]
        three_utr = [line[0], line[1], "three_prime_UTR", str(mStart), str(int(lastCDS[3]) - 1), line[5], line[6], ".", line[8]]
    print("\t".join(five_utr))
    print("\t".join(three_utr))
