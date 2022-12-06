#!/bin/bash

SN=EXAMPLE-SN
outDIR=EXAMPLE-OUT
maskDIR=/mgi_storage/sk/stomics/maskdir

docker run --rm --cpus=12 --memory=20g -u $(id -u):$(id -g) \
-v ${outDIR}:/outDir \
-v ${maskDIR}:/maskDir \
stomics/saw:04.1.0 /bin/bash CIDcount \
    -i /maskDir/${SN}.barcodeToPos.h5 \
    -s mouse \
    -g 3 \
    > ${outDIR}/00.mapping/refMemory.txt && \

# output memory estimate for mapping
# output bcNum for mapping; adde to bcPara

# genomeSize: "ls -l --block-size=GB ${Genome file after STAR indexing}"
# human/GRCh38.p13: 4
# mouse/GRCm39: 3
