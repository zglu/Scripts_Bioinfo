#!/bin/bash

inDIR=/mgi_storage/sk/stomics
outDIR=EXAMPLE-OUT
SN=SS200000183TR_D4
FC1=E100043446
FC2=E100042432
FC3=V300109183
FC4=V300109651
visualSif=/mgi_storage/sk/stomics/SAW_v5.1.3.sif

export SINGULARITY_BIND=$inDIR,$outDIR

singularity exec ${visualSif} saturation \
    -i ${outDIR}/02.count/raw_barcode_gene_exp.txt \
    --tissue ${outDIR}/04.tissuecut/${SN}.tissue.gef \
    -o ${outDIR}/06.saturation \
    --bcstat ${outDIR}/${FC1}/00.mapping/${FC1}_barcodeMap.stat,${outDIR}/${FC2}/00.mapping/${FC2}_barcodeMap.stat,${outDIR}/${FC3}/00.mapping/${FC3}_L01_barcodeMap.stat,${outDIR}/${FC3}/00.mapping/${FC3}_L02_barcodeMap.stat,${outDIR}/${FC3}/00.mapping/${FC3}_L03_barcodeMap.stat,${outDIR}/${FC3}/00.mapping/${FC3}_L04_barcodeMap.stat,${outDIR}/${FC4}/00.mapping/${FC4}_L01_barcodeMap.stat,${outDIR}/${FC4}/00.mapping/${FC4}_L02_barcodeMap.stat,${outDIR}/${FC4}/00.mapping/${FC4}_L03_barcodeMap.stat,${outDIR}/${FC4}/00.mapping/${FC4}_L04_barcodeMap.stat \
    --summary ${outDIR}/02.count/${SN}.Aligned.sortedByCoord.out.merge.q10.dedup.target.bam.summary.stat

