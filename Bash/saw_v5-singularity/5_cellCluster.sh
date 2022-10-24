#!/bin/bash

inDIR=/mgi_storage/sk/stomics
outDIR=EXAMPLE-OUT
SN=EXAMPLE-SN
visualSif=/mgi_storage/sk/stomics/SAW_v5.1.3.sif

export SINGULARITY_BIND=$inDIR,$outDIR

mkdir -p ${outDIR}/05.cellcluster

singularity exec ${visualSif} cellCluster \
    -i ${outDIR}/04.cellCut/${SN}.cellbin.gef \
    -o ${outDIR}/05.cellcluster/${SN}.cell.cluster.h5ad 
