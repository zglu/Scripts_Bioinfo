#!/bin/bash

inDIR=/mgi_storage/sk/stomics
outDIR=EXAMPLE-OUT
SN=EXAMPLE-SN
visualSif=/mgi_storage/sk/stomics/SAW_v5.4.0.sif
imgDIR=/mgi_storage/sk/stomics/ImageQC/${SN}

export SINGULARITY_BIND=$inDIR,$outDIR,$imgDIR

image4register=$(find ${imgDIR} -maxdepth 1 -name ${SN}*.tar.gz | head -1)
imageIPR=$(find ${outDIR}/03.register -maxdepth 1 -name ${SN}*.ipr | head -1)

singularity exec ${visualSif} imageTools ipr2img \
    -i ${image4register} \
    -c ${imageIPR} \
    -d tissue \
    -r True \
    -o ${outDIR}/03.register

# -d tissue cell \ if register instead of rapidRegister was done
