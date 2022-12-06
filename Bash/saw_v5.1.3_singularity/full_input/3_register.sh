#!/bin/bash


visualSif=/mgi_storage/sk/stomics/SAW_v5.1.3.sif
SN=EXAMPLE-SN
outDIR=EXAMPLE-OUT
imgDIR=/mgi_storage/sk/stomics/ssDNA/${SN}

export SINGULARITY_BIND=$outDIR,$imgDIR

imageQC=$(find ${imgDIR} -maxdepth 1 -name ${SN}*.json | head -1)
image4register=$(find ${imgDIR} -maxdepth 1 -name *.tar.gz | head -1)

singularity exec ${visualSif} register \
    -i ${image4register} \
    -c ${imageQC} \
    -v ${outDIR}/02.count/${SN}.raw.gef \
    -o ${outDIR}/03.register
