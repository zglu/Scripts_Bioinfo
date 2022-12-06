#!/bin/bash

SN=EXAMPLE-SN
outDIR=EXAMPLE-OUT

# convert tissue bin1 gef to gem
docker run --rm --cpus=64 --memory=80g -u $(id -u):$(id -g) \
-v ${outDIR}:/outDir \
stomics/saw:04.1.0 /bin/bash gefTools view \
    -i /outDir/04.tissuecut/${SN}.tissue.gef \
    -o /outDir/${SN}.tissue.gem \
    -b 1

gzip -f ${outDIR}/${SN}.tissue.gem

