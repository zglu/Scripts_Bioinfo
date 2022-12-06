#!/bin/bash

visualSif=/mgi_storage/sk/stomics/SAW_v5.1.3.sif
SN=EXAMPLE-SN
FC=EXAMPLE-FC
outDIR=EXAMPLE-OUT

export SINGULARITY_BIND=$outDIR

singularity exec ${visualSif} merge \
    --in ${outDIR}/00.mapping/${FC}_L01.barcodeReadsCount.txt,${outDIR}/00.mapping/${FC}_L02.barcodeReadsCount.txt,${outDIR}/00.mapping/${FC}_L03.barcodeReadsCount.txt,${outDIR}/00.mapping/${FC}_L04.barcodeReadsCount.txt\
    --out ${outDIR}/01.merge/${SN}.barcodeReadsCount.txt \
    --action 2

