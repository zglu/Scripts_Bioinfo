#!/bin/bash

inDIR=/mgi_storage/sk/stomics
SN=EXAMPLE-SN
FC=EXAMPLE-FC
outDIR=EXAMPLE-OUT
visualSif=/mgi_storage/sk/stomics/SAW_v5.1.3.sif

export SINGULARITY_BIND=$inDIR,$outDIR

singularity exec ${visualSif} merge \
    ${inDIR}/maskdir/${SN}.barcodeToPos.h5\
    ${outDIR}/00.mapping/${FC}_L01.barcodeReadsCount.txt,${outDIR}/00.mapping/${FC}_L02.barcodeReadsCount.txt,${outDIR}/00.mapping/${FC}_L03.barcodeReadsCount.txt,${outDIR}/00.mapping/${FC}_L04.barcodeReadsCount.txt\
    ${outDIR}/01.merge/${SN}.barcodeReadsCount.txt
