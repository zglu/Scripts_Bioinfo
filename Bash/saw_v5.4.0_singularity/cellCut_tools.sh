#!/bin/bash

inDIR=/mgi_storage/sk/stomics
outDIR=EXAMPLE-OUT
SN=EXAMPLE-SN
visualSif=/mgi_storage/sk/stomics/SAW_v5.4.0.sif

export SINGULARITY_BIND=$inDIR,$outDIR

# convert tissue bin1 gef to gem
singularity exec ${visualSif} cellCut view \
    -s ${SN} \
    -i ${outDIR}/04.tissuecut/${SN}.tissue.gef \
    -o ${outDIR}/analysis/${SN}.tissue.gem \

# convert whole chip bin1 gef to gem
singularity exec ${visualSif} cellCut view \
    -s ${SN} \
    -i ${outDIR}/02.count/${SN}.raw.gef \
    -o ${outDIR}/analysis/${SN}.raw.gem 

# convert cellbin bin1 gef to bem
#singularity exec ${visualSif} cellCut view \
#    -s ${SN} \
#    -i ${outDIR}/04.cellcut/${SN}.cellbin.gef \
#    -o ${outDIR}/analysis/${SN}.cellbin.gem \
#    -d ${outDIR}/02.count/${SN}.raw.gef

#gzip -f ${outDIR}/analysis/${SN}.tissue.bin1.gem

