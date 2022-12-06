#!/bin/bash

visualSif=/mgi_storage/sk/stomics/SAW_v5.1.3.sif
SN=EXAMPLE-SN
FC=EXAMPLE-FC # E100010821_L01_120
outDIR=EXAMPLE-OUT
fqDIR=EXAMPLE-FQ #
maskDIR=/mgi_storage/sk/stomics/maskdir
refDIR=/mgi_storage/sk/stomics/reference/GRCm39_GENCODE

export SINGULARITY_BIND=$outDIR,$fqDIR,$maskDIR,$refDIR

singularity exec ${visualSif} mapping \
    --outSAMattributes spatial \
    --outSAMtype BAM SortedByCoordinate \
    --genomeDir ${refDIR}/STAR_SJ100 \
    --runThreadN 36 \
    --outFileNamePrefix ${outDIR}/00.mapping/${FC}_${LANE}. \
    --sysShell /bin/bash \
    --stParaFile ${outDIR}/00.mapping/L01.bcPara \
    --readNameSeparator \" \" \
    --limitBAMsortRAM 80000000000 \
    --limitOutSJcollapsed 10000000 \
    --limitIObufferSize=280000000 \
    --outBAMsortingBinsN 50 \
    > ${outDIR}/00.mapping/${FC}_${LANE}_barcodeMap.stat

# human: refDIR=/mgi_storage/sk/stomics/reference/GRCh38.p13_GENCODE
# mouse: refDIR=/mgi_storage/sk/stomics/reference/GRCm39_GENCODE
# wheet: refDIR=/mgi_storage/sk/stomics/reference/Triticum_aestivum.IWGSC54 Triticum_aestivum.IWGSC.54.gff3
# EXITING: FATAL INPUT ERROR 可从命令中删除 --readNameSeparator \" \":
# multiple lanes: need to run mapping for each lane
