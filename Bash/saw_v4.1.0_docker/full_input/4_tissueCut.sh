#!/bin/bash

SN=EXAMPLE-SN
outDIR=EXAMPLE-OUT

docker run --rm --cpus=24 --memory=60g -u $(id -u):$(id -g) \
  -v ${outDIR}:/outDir \
stomics/saw:04.1.0 /bin/bash tissueCut \
    --dnbfile /outDir/01.merge/${SN}.barcodeReadsCount.txt \
    -i /outDir/02.count/${SN}.raw.gef \
    -o /outDir/04.tissuecut \
#   -s /outDir/04.register/7_result \
    -t tissue \
    --platform T10 \
    --snId ${SN}

