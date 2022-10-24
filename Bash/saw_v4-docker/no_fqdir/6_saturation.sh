#!/bin/bash

SN=EXAMPLE-SN
FC=EXAMPLE-FC
outDIR=EXAMPLE-OUT

docker run --rm --cpus=24 --memory=60g \
  -v ${outDIR}:/outDir \
stomics/saw:04.1.0 /bin/bash saturation \
    -i /outDir/02.count/raw_barcode_gene_exp.txt \
    --tissue /outDir/04.tissuecut/${SN}.tissue.gef \
    -o /outDir/06.saturation \
    --bcstat /outDir/00.mapping/${FC}_L01_barcodeMap.stat,/outDir/00.mapping/${FC}_L02_barcodeMap.stat,/outDir/00.mapping/${FC}_L03_barcodeMap.stat,/outDir/00.mapping/${FC}_L04_barcodeMap.stat \
    --summary /outDir/02.count/${SN}.Aligned.sortedByCoord.out.merge.q10.dedup.target.bam.summary.stat &&\

