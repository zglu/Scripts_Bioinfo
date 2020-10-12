#! /usr/bin/env python

"""
Convert a gff3 annotation file to a sqlite3 database file
zl3 // last update: 11/10/2017
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
myDB = 'annotation.db'

db = gffutils.create_db(myGFF, myDB, merge_strategy="create_unique", keep_order=True) ## :memory:

## or create a real database
## import data into database: one-time operation
#db = gffutils.create_db(myGFF, 'myGFF.db', merge_strategy="create_unique", keep_order=True) # ':memory:'
## once created, use the database (donot need gff file anymore)
#db = gffutils.FeatureDB('myGFF.db')

