#! /usr/bin/env python

"""
Get dna sequences based on a gff file for features like gene, transcript, cds, exon, and get pep for coding sequences.
SHOULD CHECK THE FEATURE NAME IN YOUR GFF: mRNA or transcript, and change accordingly.
Usage: python gff2fasta.py <gff3 file> <genome sequence> <feature> > output
last update: 17/07/2017
"""

import sys
import gffutils
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna

## ignore biopython warnings
import warnings
from Bio import BiopythonWarning
warnings.simplefilter('ignore', BiopythonWarning)

myGFF = sys.argv[1]
myFasta = sys.argv[2]

db = gffutils.create_db(myGFF, ':memory:', merge_strategy="create_unique", keep_order=True)

## function to get all gene/transcript(mRNA) sequences
def parent_seq(type):
    for p in db.features_of_type(type):
        p_seq = p.sequence(myFasta)
        p_seq = Seq(p_seq, generic_dna)
        if p.strand == '-':
            p_seq = p_seq.reverse_complement()
        print('>' + p.id + '\n' + p_seq)

## function to get all cds/exon/intron sequences under the same transcript id
def child_seq(type):
    for t in db.features_of_type('transcript', order_by='start'): # or mRNA depending on the gff
        # print(t.sequence(myFasta))
        seq_combined = ''
        for i in db.children(t, featuretype=type, order_by='start'): # or exon/intron
            seq = i.sequence(myFasta, use_strand=True)  # use_strand doesn't work; have to revcomp
            seq_combined += seq
        seq_combined = Seq(seq_combined, generic_dna)
        if t.strand == '-':
            seq_combined = seq_combined.reverse_complement()
        print('>' + t.id + '\n' + seq_combined)

## function to get all protein sequences under the same transcript id
def pep_seq():
    for t in db.features_of_type('transcript', order_by='start'): # or mRNA depending on the gff
        # print(t.sequence(myFasta))
        seq_combined = ''
        for i in db.children(t, featuretype='CDS', order_by='start'): # or exon/intron
            seq = i.sequence(myFasta, use_strand=True)  # use_strand doesn't work; have to revcomp
            seq_combined += seq
        seq_combined = Seq(seq_combined, generic_dna)
        if t.strand == '-':
            seq_combined = seq_combined.reverse_complement()
        print('>' + t.id + '\n' + seq_combined.translate())

seqType = sys.argv[3]

if seqType in ['gene', 'transcript']:
    parent_seq(seqType)
elif seqType in ['CDS', 'exon', 'intron']:
    child_seq(seqType)
elif seqType == 'pep':
    pep_seq()
else:
    sys.exit("Supported feature types: gene/transcript/CDS/pep/exon/intron.")