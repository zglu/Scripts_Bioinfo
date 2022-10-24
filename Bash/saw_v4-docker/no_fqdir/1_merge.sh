#!/bin/bash

SN=EXAMPLE-SN
FC=EXAMPLE-FC
outDIR=EXAMPLE-OUT

docker run --rm --cpus=12 --memory=20g \
  -v ${outDIR}:/outDir \
stomics/saw:04.1.0 /bin/bash merge \
    --in /outDir/00.mapping/${FC}_L01.barcodeReadsCount.txt,/outDir/00.mapping/${FC}_L02.barcodeReadsCount.txt,/outDir/00.mapping/${FC}_L03.barcodeReadsCount.txt,/outDir/00.mapping/${FC}_L04.barcodeReadsCount.txt\
    --out /outDir/01.merge/${SN}.barcodeReadsCount.txt \
    --action 2 && \

