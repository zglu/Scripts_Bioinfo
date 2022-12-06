#!/bin/bash

inDIR=/mgi_storage/sk/stomics
outDIR=EXAMPLE-OUT
SN=EXAMPLE-SN
visualSif=/mgi_storage/sk/stomics/SAW_v5.1.3.sif

export SINGULARITY_BIND=$inDIR,$outDIR

singularity exec ${visualSif} tissueCut \
    --dnbfile ${outDIR}/01.merge/${SN}.barcodeReadsCount.txt \
    -i ${outDIR}/02.count/${SN}.raw.gef \
    -o ${outDIR}/04.tissuecut \
    -s ${outDIR}/03.register \
    -t tissue \
    --snId ${SN} \
    --omics=Transcriptomics \
    -d

# delete -s... line if no register was done.
