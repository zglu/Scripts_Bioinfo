#!/bin/bash

SN=EXAMPLE-SN
outDIR=EXAMPLE-OUT

docker run --rm --cpus=12 --memory=20g \
  -v ${outDIR}:/outDir \
stomics/saw:04.1.0 /bin/bash spatialCluster \
    -i /outDir/04.tissuecut/${SN}.tissue.gef \
    -o /outDir/05.spatialcluster/${SN}.spatial.cluster.h5ad \
    -s 50 &&\

# -s: binSize
