#!/bin/bash

inDIR=/mgi_storage/sk/stomics
inDIR2=/mgi_storage/customer/mgi_lvprod/share/mosta/E16.5_E1S3
SN=EXAMPLE-SN
FC=EXAMPLE-FC 
SB=EXAMPLE-SB
outDIR=EXAMPLE-OUT
visualSif=/mgi_storage/sk/stomics/SAW_v5.4.0.sif

export SINGULARITY_BIND=$inDIR,$inDIR2,$outDIR

singularity exec ${visualSif} mapping \
    --outSAMattributes spatial \
    --outSAMtype BAM SortedByCoordinate \
    --genomeDir ${inDIR}/reference/GRCm39_GENCODE/STAR_SJ100 \
    --runThreadN 36 \
    --outFileNamePrefix ${outDIR}/00.mapping/${FC}_${SB}. \
    --sysShell /bin/bash \
    --stParaFile ${outDIR}/00.mapping/L01_${SB}.bcPara \
    --readNameSeparator \" \" \
    --limitBAMsortRAM 200000000000 \
    --limitOutSJcollapsed 10000000 \
    --limitIObufferSize=280000000 \
    --outBAMsortingBinsN 50 \
    > ${outDIR}/00.mapping/${FC}_${SB}_barcodeMap.stat

# EXITING: FATAL INPUT ERROR 可从命令中删除 --readNameSeparator \" \":
