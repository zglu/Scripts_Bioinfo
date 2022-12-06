#!/bin/bash

inDIR=/mgi_storage/sk/stomics
outDIR=EXAMPLE-OUT
SN=EXAMPLE-SN
visualSif=/mgi_storage/sk/stomics/SAW_v5.4.0.sif

export SINGULARITY_BIND=$inDIR,$outDIR

tissueMaskFile=$(find ${outDIR}/03.register -maxdepth 1 -name {SN}_tissue_cut.tif | head -1)

singularity exec ${visualSif} tissueCut \
    --dnbfile ${outDIR}/01.merge/${SN}.barcodeReadsCount.txt \
    -i ${outDIR}/02.count/${SN}.raw.gef \
    -o ${outDIR}/04.tissuecut \
    -s ${tissueMaskFile} \
    --sn ${SN} \
    --omics=Transcriptomics \
    -d

cp stereo-chip.tissue.gef ${SN}.tissue.gef

# delete -s... if no register was done
