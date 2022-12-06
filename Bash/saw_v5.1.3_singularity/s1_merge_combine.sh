#!/bin/bash

inDIR=/mgi_storage/sk/stomics
SN=SS200000183TR_D4
FC1=E100043446
FC2=E100042432
FC3=V300109183
FC4=V300109651
outDIR=EXAMPLE-OUT
visualSif=/mgi_storage/sk/stomics/SAW_v5.1.3.sif

export SINGULARITY_BIND=$inDIR,$outDIR

singularity exec ${visualSif} merge \
    ${inDIR}/maskdir/${SN}.barcodeToPos.h5\
    ${outDIR}/${FC1}/00.mapping/${FC1}.barcodeReadsCount.txt,${outDIR}/${FC2}/00.mapping/${FC2}.barcodeReadsCount.txt,${outDIR}/${FC3}/00.mapping/${FC3}_L01.barcodeReadsCount.txt,${outDIR}/${FC3}/00.mapping/${FC3}_L02.barcodeReadsCount.txt,${outDIR}/${FC3}/00.mapping/${FC3}_L03.barcodeReadsCount.txt,${outDIR}/${FC3}/00.mapping/${FC3}_L04.barcodeReadsCount.txt,${outDIR}/${FC4}/00.mapping/${FC4}_L01.barcodeReadsCount.txt,${outDIR}/${FC4}/00.mapping/${FC4}_L02.barcodeReadsCount.txt,${outDIR}/${FC4}/00.mapping/${FC4}_L03.barcodeReadsCount.txt,${outDIR}/${FC4}/00.mapping/${FC4}_L04.barcodeReadsCount.txt \
    ${outDIR}/01.merge/${SN}.barcodeReadsCount.txt

