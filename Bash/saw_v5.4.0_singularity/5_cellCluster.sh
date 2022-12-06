#!/bin/bash

inDIR=/mgi_storage/sk/stomics
outDIR=EXAMPLE-OUT
SN=EXAMPLE-SN
visualSif=/mgi_storage/sk/stomics/SAW_v5.4.0.sif

export SINGULARITY_BIND=$inDIR,$outDIR

mkdir -p ${outDIR}/05.cellcluster

singularity exec ${visualSif} cellCluster \
    -i ${outDIR}/04.cellcut/${SN}.cellbin.gef \
    -o ${outDIR}/05.cellcluster/${SN}.cell.cluster.h5ad
