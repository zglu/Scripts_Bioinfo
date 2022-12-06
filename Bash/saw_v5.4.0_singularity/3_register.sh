#!/bin/bash

inDIR=/mgi_storage/sk/stomics
outDIR=EXAMPLE-OUT
SN=EXAMPLE-SN
visualSif=/mgi_storage/sk/stomics/SAW_v5.4.0.sif
imgDIR=/mgi_storage/sk/stomics/ImageQC/${SN}

export SINGULARITY_BIND=$inDIR,$outDIR,$imgDIR

imageQC=$(find ${imgDIR} -maxdepth 1 -name ${SN}*.ipr | head -1)
image4register=$(find ${imgDIR} -maxdepth 1 -name *.tar.gz | head -1)

singularity exec ${visualSif} register \
    -i ${image4register} \
    -c ${imageQC} \
    -v ${outDIR}/02.count/${SN}.raw.gef \
    -o ${outDIR}/03.register
