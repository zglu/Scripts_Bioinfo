#!/bin/bash

visualSif=/mgi_storage/sk/stomics/SAW_v5.1.3.sif
SN=EXAMPLE-SN
outDIR=EXAMPLE-OUT
maskDIR=/mgi_storage/sk/stomics/maskdir

export SINGULARITY_BIND=$outDIR,$maskDIR

singularity exec ${visualSif} CIDcount \
    -i ${maskDIR}/${SN}.barcodeToPos.h5 \
    -s mouse \
    -g 3 \
    > ${outDIR}/00.mapping/refMemory.txt

# output memory estimate for mapping
# output bcNum for mapping; adde to bcPara

# genomeSize: "ls -l --block-size=GB ${Genome file after STAR indexing}"
# human/GRCh38.p13: 4
# mouse/GRCm39: 3
