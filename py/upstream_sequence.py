#! /usr/bin/env python

"""
Get upstream sequences of gene/mRNA features from one gff file.
last update: 28/06/2018
"""

import sys
import gffutils
from Bio.Seq import Seq
from pyfasta import Fasta

## ignore biopython warnings
import warnings
from Bio import BiopythonWarning
warnings.simplefilter('ignore', BiopythonWarning)

if len(sys.argv) < 4:
    sys.exit("Usage: python upstream_seq.py <gff> <fasta> <feature> <range> > output")
else:
    myGFF = sys.argv[1]
    myFasta = sys.argv[2]
    myType = sys.argv[3]
    myRange = int(sys.argv[4])

db = gffutils.create_db(myGFF, ':memory:', merge_strategy="create_unique", keep_order=True)

def upstream_seq(type):
    for p in db.features_of_type(type):
        genomefa = Fasta(myFasta)
        print('>' + p.id + "_Upstream" + str(myRange))
        upstart = p.start - 1 - myRange
        upend = p.start -1
        # get sequence based on coordinates (start is 0-based)
        p_upstream = genomefa[p.seqid][upstart:upend]
        if p.strand == '-':
            upstart = p.end
            upend = p.end + myRange
            p_upstream = genomefa[p.seqid][upstart:upend]
            p_upstream = Seq(p_upstream).reverse_complement()
        for i in range(0, len(p_upstream), 60): # print 60 bases per line
            print(p_upstream[i:i+60])

if myType in ['gene', 'mRNA']:
    upstream_seq(myType)
else:
    sys.exit("Feature should be gene or mRNA")