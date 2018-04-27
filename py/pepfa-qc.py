#! /usr/bin/env python

"""
QC for a peptide fasta file. Will report the Start, End and number of internal stop codons for each sequence.
"""

import sys
import re
from Bio import SeqIO
#import numpy as np
#import matplotlib.pyplot as plt

myFasta = sys.argv[1]

fa_sequences = SeqIO.parse(myFasta, 'fasta')
#shortseq_list = [fa for fa in fa_sequences if len(fa.seq) < 50]
#shortseq_generator = (fa for fa in fa_sequences if len(fa.seq) < 50)
#SeqIO.write(shortseq_generator, "short_sequences.fasta", "fasta")

with open("pepfa-qc.txt", 'w') as fout:
    for fa in fa_sequences:
        id, sequence = fa.id, str(fa.seq)
        startc = sequence[0]
        endc = sequence[-1]
        errs = re.findall('\*', sequence[:-1])
        fout.write(id + '\tLEN= ' + str(len(sequence)-1) + '\tSTART= ' + startc + '\tEND= ' + endc + '\tINTERNAL_STOPS= ' + str(len(errs)) + '\n')
fout.close()
