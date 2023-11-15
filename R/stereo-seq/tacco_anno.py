
import os
import sys
import matplotlib

import pandas as pd
import numpy as np
import anndata as ad

import tacco as tc

scref=sys.argv[1]
spdata=sys.argv[2]
grpname=sys.argv[3]

reference = ad.read(scref)
puck=ad.read(spdata)

newKey='tacco_' + grpname

tc.tl.annotate(puck, reference, grpname, result_key=newKey)
puck.obsm[newKey].to_csv(spdata + "_" + newKey + ".csv")


#fig = tc.pl.scatter(puck, keys=newKey, position_key=['coord_x','coord_y'], joint=True, point_size=1, 
#                    noticks=True, axes_labels=['X','Y']);
