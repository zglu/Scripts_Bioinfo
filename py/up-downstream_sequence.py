#! /usr/bin/env python

"""
Get up- and downstream (relative to transcription start site) sequences of gene/mRNA features from one gff file.
last update: 12/06/2019
"""

import sys
import gffutils
from Bio.Seq import Seq
from pyfasta import Fasta

## ignore biopython warnings
import warnings
from Bio import BiopythonWarning
warnings.simplefilter('ignore', BiopythonWarning)

if len(sys.argv) < 5:
    sys.exit("Usage: python3 up-downstream_seq.py <gff> <fasta> <feature> <upRange> <downRange> > output")
else:
    myGFF = sys.argv[1]
    myFasta = sys.argv[2]
    myType = sys.argv[3]
    upRange = int(sys.argv[4])
    downRange = int(sys.argv[5])

db = gffutils.create_db(myGFF, ':memory:', merge_strategy="create_unique", keep_order=True)

def updownTSS_seq(type):
    for p in db.features_of_type(type):
        genomefa = Fasta(myFasta)
        print('>' + p.id + "_[-" + str(upRange) + "]-[+" + str(downRange) + "]")
        seqstart = p.start - 1 - upRange
        seqstart = seqstart if seqstart > 0 else 1  # avoid start with minus coord
        seqend = p.start + downRange
        # get sequence based on coordinates (start is 0-based)
        p_updown = genomefa[p.seqid][seqstart:seqend]
        if p.strand == '-':
            seqstart = p.end - 1 - downRange
            seqstart = seqstart if seqstart > 0 else 1 # avoid start with minus coord 
            seqend = p.end + upRange
            p_updown = genomefa[p.seqid][seqstart:seqend]
            p_updown = Seq(p_updown).reverse_complement()
        for i in range(0, len(p_updown), 60): # print 60 bases per line
            print(p_updown[i:i+60])

if myType in ['gene', 'mRNA']:
    updownTSS_seq(myType)
else:
    sys.exit("Feature should be gene or mRNA")
