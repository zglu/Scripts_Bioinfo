#! /usr/bin/env python

import sys
from Bio.SeqIO.FastaIO import SimpleFastaParser

FastaFile = open(sys.argv[1], 'rU')

for name, seq in SimpleFastaParser(FastaFile):
    seqLen = len(seq)
    print(name + '\t' + str(seqLen))

FastaFile.close()
