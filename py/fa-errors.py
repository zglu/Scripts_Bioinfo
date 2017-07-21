#! /usr/bin/env python

"""
Find number of stop codons in the protein sequences.
"""

import sys
import re
from Bio import SeqIO
import numpy as np
import matplotlib.pyplot as plt

myFasta = sys.argv[1]

fa_sequences = SeqIO.parse(myFasta, 'fasta')
#shortseq_list = [fa for fa in fa_sequences if len(fa.seq) < 50]
#shortseq_generator = (fa for fa in fa_sequences if len(fa.seq) < 50)
#SeqIO.write(shortseq_generator, "short_sequences.fasta", "fasta")

#with open("errlist.txt", 'w') as fout:
#    for fa in fa_sequences:
#        id, sequence = fa.id, str(fa.seq)
#        errlist = re.findall('\*', sequence)
#        if len(errlist) > 2:
#            fout.write(id + ' ' + str(len(errlist)) + '\n')
#fout.close()

errList = []
for fa in fa_sequences:
    id, sequence = fa.id, str(fa.seq)
    errs = re.findall('\*', sequence)
    errList.append(len(errs))
#    if len(errlist) > 20:
#        print(id + ': ' + str(len(errs)) + ' Stop Codons in the sequence.')
#print(errList)
print(len(errList))
errArray = np.asarray(errList)
#newArray = np.array(errArray).astype(np.int)

## https://docs.scipy.org/doc/numpy/reference/routines.statistics.html
print('min: ')
print(np.amin(errArray))
print('max: ')
print(np.amax(errArray))

hist, bin_edges = np.histogram(errArray, bins=[1, 10, 50, 100, 120]) # intervals
print(bin_edges)
print(hist)

## show the histogram
#plt.hist(errArray, bins=[0, 10, 20, 50, 100, 120])
#plt.show()