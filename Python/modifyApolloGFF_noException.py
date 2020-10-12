#! /usr/bin/env python

'''
modify GFF from Apollo to Augustus style
must be in gene -> mRNA -> cds/exon order
first remove comment and empty lines: grep -v '#' | grep -v -e '^$'
'''

# last update: 21/12/2017

import re
import sys

GFF = sys.argv[1] # gff with many genes


with open(GFF) as fin:
    for line in fin:
        line = line.rstrip()
        line = line.split()
        if "gene" in line or "pseudogene" in line:
            i = 0
            name_pat = re.compile("Name=(Smp_......);")
            geneid = name_pat.search(line[8]).group(1)
            line[8] = "ID=" + geneid
        elif "mRNA" in line or "transcript" in line:
            mrnaid_pat = re.compile("ID=(.+?);")
            mrnaid = mrnaid_pat.search(line[8]).group(1)
            i += 1
            trid = geneid + "." + str(i)
            line[8] = "ID=" + trid + ";" + "Parent=" + geneid
        elif any(s in line for s in ("exon", "non_canonical_three_prime_splice_site", "non_canonical_five_prime_splice_site")):
            line[8] = line[8].replace(mrnaid, trid)
            line[8] = re.sub(";ID=(.+)", "", line[8])
            line[8] = re.sub("Name=(.+)", "", line[8])
        elif "CDS" in line:
            line[8] = line[8].replace(mrnaid, trid)
            line[8] = re.sub("Parent=(.+?);", "", line[8])
            line[8] = re.sub("ID=(.+);", "ID=" + trid + ".cds;Parent=" + trid, line[8])
            line[8] = re.sub("Name=(.+)", "", line[8])
        line[1] = "Apollo"
        print("\t".join(line))
