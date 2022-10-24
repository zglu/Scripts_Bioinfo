#!/bin/bash

inDIR=/mgi_storage/sk/stomics
outDIR=EXAMPLE-OUT
SN=EXAMPLE-SN
visualSif=/mgi_storage/sk/stomics/SAW_v5.1.3.sif

export SINGULARITY_BIND=$inDIR,$outDIR

singularity exec ${visualSif} CIDCount \
    -i ${inDIR}/maskdir/${SN}.barcodeToPos.h5 \
    -s mouse \
    -g 3 \
    > ${outDIR}/00.mapping/refMemory.txt 

# output memory estimate for mapping
# output bcNum for mapping; adde to bcPara

# genomeSize: "ls -l --block-size=GB ${Genome file after STAR indexing}"
# human/GRCh38.p13: 4
# mouse/GRCm39: 3
