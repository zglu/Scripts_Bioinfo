#!/bin/bash

inDIR=/mgi_storage/sk/stomics
outDIR=EXAMPLE-OUT
SN=EXAMPLE-SN
visualSif=/mgi_storage/sk/stomics/SAW_v5.1.3.sif

export SINGULARITY_BIND=$inDIR,$outDIR

singularity exec ${visualSif} spatialCluster \
    -i ${outDIR}/04.tissuecut/${SN}.tissue.gef \
    -o ${outDIR}/05.spatialcluster/${SN}.spatial.cluster.h5ad \
    -s 200

# -s: binSize
