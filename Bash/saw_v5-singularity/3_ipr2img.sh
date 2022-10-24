#!/bin/bash

inDIR=/mgi_storage/sk/stomics
outDIR=EXAMPLE-OUT
SN=EXAMPLE-SN
visualSif=/mgi_storage/sk/stomics/SAW_v5.1.3.sif
imgDIR=/mgi_storage/sk/stomics/ssDNA/${SN}

export SINGULARITY_BIND=$inDIR,$outDIR,$imgDIR

image4register=$(find ${imgDIR} -maxdepth 1 -name ${SN}*.tar.gz | head -1)
imageIPR=$(find ${outDIR}/03.register -maxdepth 1 -name ${SN}*.ipr | head -1)

singularity exec ${visualSif} ipr2img \
    -i ${image4register} \
    -c ${imageIPR} \
    -d tissue cell \
    -r True \
    -o ${outDIR}/03.register
