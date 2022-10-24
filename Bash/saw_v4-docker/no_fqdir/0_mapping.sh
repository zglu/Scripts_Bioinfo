#!/bin/bash

# --SHOULD ADD bcPara FILES BEFORE MAPPING--

SN=EXAMPLE-SN
FC=EXAMPLE-FC
LANE=L01
outDIR=EXAMPLE-OUT
fqDIR=/mgi_storage/sk/stomics/fastq/${FC}
maskDIR=/mgi_storage/sk/stomics/maskdir
refDIR=/mgi_storage/sk/stomics/reference/GRCm39_GENCODE

docker run --rm --cpus=64 --memory=80g \
  -v ${outDIR}:/outDir \
  -v ${fqDIR}:/fqDir \
  -v ${maskDIR}:/maskDir \
  -v ${refDIR}:/refDir \
stomics/saw:04.1.0 /bin/bash mapping \
    --outSAMattributes spatial \
    --outSAMtype BAM SortedByCoordinate \
    --genomeDir /refDir/STAR_SJ100 \
    --runThreadN 36 \
    --outFileNamePrefix /outDir/00.maping/${FC}_${LANE}. \
    --sysShell /bin/bash \
    --stParaFile /outDir/00.mapping/${LANE}.bcPara \
    --readNameSeparator \" \" \
    --limitBAMsortRAM 63168332971 \
    --limitOutSJcollapsed 10000000 \
    --limitIObufferSize=280000000 \
    --outBAMsortingBinsN 50 \
    > ${outDIR}/00.mapping/${FC}_${LANE}_barcodeMap.stat &&\

# human: refDIR=/mgi_storage/sk/stomics/reference/GRCh38.p13_GENCODE
# mouse: refDIR=/mgi_storage/sk/stomics/reference/GRCm39_GENCODE
# EXITING: FATAL INPUT ERROR 可从命令中删除 --readNameSeparator \" \":
# multiple lanes: need to run mapping for each lane
