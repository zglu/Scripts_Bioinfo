#! /usr/bin/env python3

"""
Split a fasta file with x sequences per file
From http://biopython.org/wiki/Split_large_file; need to change iterator.next() to next(iterator)
"""

import sys
from Bio import SeqIO

if len(sys.argv) < 2:
    sys.exit("Usage: python fa-split.py FASTA_FILE SPLIT_SIZE")
else:
    myFasta = sys.argv[1]
    split_size = sys.argv[2] # e.g. 1000 sequences per file

def batch_iterator(iterator, batch_size):
    """Returns lists of length batch_size.

    This can be used on any iterator, for example to batch up
    SeqRecord objects from Bio.SeqIO.parse(...), or to batch
    Alignment objects from Bio.AlignIO.parse(...), or simply
    lines from a file handle.

    This is a generator function, and it returns lists of the
    entries from the supplied iterator.  Each list will have
    batch_size entries, although the final list may be shorter.
    """
    entry = True  # Make sure we loop once
    while entry:
        batch = []
        while len(batch) < batch_size:
            try:
                entry = next(iterator)
            except StopIteration:
                entry = None
            if entry is None:
                # End of file
                break
            batch.append(entry)
        if batch:
            yield batch

record_iter = SeqIO.parse(open(myFasta),"fasta")
for i, batch in enumerate(batch_iterator(record_iter, split_size)):
    filename = "Part-%i.fasta" % (i + 1)
    with open(filename, "w") as handle:
        count = SeqIO.write(batch, handle, "fasta")
    print("Wrote %i records to %s" % (count, filename))