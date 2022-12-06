#!/bin/bash

visualSif=/mgi_storage/sk/stomics/SAW_v5.1.3.sif
SN=EXAMPLE-SN
outDIR=EXAMPLE-OUT

export SINGULARITY_BIND=$outDIR

# convert tissue bin1 gef to gem
singularity exec ${visualSif} gefTools view \
    -i ${outDIR}/04.tissuecut/${SN}.tissue.gef \
    -o ${outDIR}/${SN}.tissue.gem \
    -b 1

gzip -f ${outDIR}/${SN}.tissue.gem

