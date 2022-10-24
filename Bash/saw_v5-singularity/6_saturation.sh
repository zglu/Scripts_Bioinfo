#!/bin/bash

inDIR=/mgi_storage/sk/stomics
outDIR=EXAMPLE-OUT
SN=EXAMPLE-SN
FC=EXAMPLE-FC
visualSif=/mgi_storage/sk/stomics/SAW_v5.1.3.sif

export SINGULARITY_BIND=$inDIR,$outDIR

singularity exec ${visualSif} saturation \
    -i ${outDIR}/02.count/raw_barcode_gene_exp.txt \
    --tissue ${outDIR}/04.tissuecut/${SN}.tissue.gef \
    -o ${outDIR}/06.saturation \
    --bcstat ${outDIR}/00.mapping/${FC}_L01_barcodeMap.stat,${outDIR}/00.mapping/${FC}_L02_barcodeMap.stat,${outDIR}/00.mapping/${FC}_L03_barcodeMap.stat,${outDIR}/00.mapping/${FC}_L04_barcodeMap.stat \
    --summary ${outDIR}/02.count/${SN}.Aligned.sortedByCoord.out.merge.q10.dedup.target.bam.summary.stat

