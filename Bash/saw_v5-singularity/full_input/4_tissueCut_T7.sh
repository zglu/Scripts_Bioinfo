#!/bin/bash

visualSif=/mgi_storage/sk/stomics/SAW_v5.1.3.sif
SN=EXAMPLE-SN
FC=EXAMPLE-FC
outDIR=EXAMPLE-OUT

export SINGULARITY_BIND=$outDIR

singularity exec ${visualSif} tissueCut \
    --dnbfile ${outDIR}/00.mapping/${FC}.barcodeReadsCount.txt \
    -i ${outDIR}/02.count/${SN}.raw.gef \
    -o ${outDIR}/04.tissuecut \
#   -s ${outDIR}/04.register/7_result \
    -t tissue \
    --platform T10 \
    --snId ${SN}

