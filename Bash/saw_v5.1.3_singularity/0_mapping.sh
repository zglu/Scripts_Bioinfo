#!/bin/bash

inDIR=/mgi_storage/sk/stomics
outDIR=EXAMPLE-OUT
SN=EXAMPLE-SN
FC=EXAMPLE-FC
LANE=L01
visualSif=/mgi_storage/sk/stomics/SAW_v5.1.3.sif

export SINGULARITY_BIND=$inDIR,$outDIR

singularity exec ${visualSif} mapping \
    --outSAMattributes spatial \
    --outSAMtype BAM SortedByCoordinate \
    --genomeDir ${inDIR}/reference/GRCm39_GENCODE/STAR_SJ100 \
    --runThreadN 36 \
    --outFileNamePrefix ${outDIR}/00.mapping/${FC}_${LANE}. \
    --sysShell /bin/bash \
    --stParaFile ${outDIR}/00.mapping/${LANE}.bcPara \
    --readNameSeparator \" \" \
    --limitBAMsortRAM 80000000000 \
    --limitOutSJcollapsed 10000000 \
    --limitIObufferSize=280000000 \
    --outBAMsortingBinsN 50 \
    > ${outDIR}/00.mapping/${FC}_${LANE}_barcodeMap.stat

# EXITING: FATAL INPUT ERROR 可从命令中删除 --readNameSeparator \" \":
