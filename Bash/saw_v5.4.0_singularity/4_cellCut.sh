#!/bin/bash

inDIR=/mgi_storage/sk/stomics
outDIR=EXAMPLE-OUT
SN=EXAMPLE-SN
visualSif=/mgi_storage/sk/stomics/SAW_v5.4.0.sif

export SINGULARITY_BIND=$inDIR,$outDIR

mkdir -p ${outDIR}/04.cellcut

singularity exec ${visualSif} cellCut cgef \
    -i ${outDIR}/02.count/${SN}.raw.gef \
    -m ${outDIR}/03.register/${SN}_mask.tif \
    -o ${outDIR}/04.cellcut/${SN}.cellbin.gef
