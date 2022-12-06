#!/bin/bash

visualSif=/mgi_storage/sk/stomics/SAW_v5.1.3.sif
SN=EXAMPLE-SN
outDIR=EXAMPLE-OUT

export SINGULARITY_BIND=$outDIR

singularity exec ${visualSif} spatialCluster \
    -i ${outDIR}/04.tissuecut/${SN}.tissue.gef \
    -o ${outDIR}/05.spatialcluster/${SN}.spatial.cluster.h5ad \
    -s 50

# -s: binSize
