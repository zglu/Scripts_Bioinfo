#! /usr/bin/env python

"""
Get up- and downstream sequences of gene/mRNA features from one gff file.
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

def upstream_seq(type):
    for p in db.features_of_type(type):
        genomefa = Fasta(myFasta)
        print('>' + p.id + "_[-" + str(upRange) + "]-[+" + str(downRange) + "]")
        upstart = p.start - 1 - upRange
        upend = p.start + downRange
        # get sequence based on coordinates (start is 0-based)
        p_upstream = genomefa[p.seqid][upstart:upend]
        if p.strand == '-':
            upstart = p.end - 1 - downRange
            upend = p.end + upRange
            p_upstream = genomefa[p.seqid][upstart:upend]
            p_upstream = Seq(p_upstream).reverse_complement()
        for i in range(0, len(p_upstream), 60): # print 60 bases per line
            print(p_upstream[i:i+60])

if myType in ['gene', 'mRNA']:
    upstream_seq(myType)
else:
    sys.exit("Feature should be gene or mRNA")
