#!/bin/bash

visualSif=/mgi_storage/sk/stomics/SAW_v5.1.3.sif
SN=EXAMPLE-SN
FC=EXAMPLE-FC
outDIR=EXAMPLE-OUT

export SINGULARITY_BIND=$outDIR

singularity exec ${visualSif} saturation \
    -i ${outDIR}/02.count/raw_barcode_gene_exp.txt \
    --tissue ${outDIR}/04.tissuecut/${SN}.tissue.gef \
    -o ${outDIR}/06.saturation \
    --bcstat ${outDIR}/00.mapping/${FC}_barcodeMap.stat \
    --summary ${outDIR}/02.count/${SN}.Aligned.sortedByCoord.out.merge.q10.dedup.target.bam.summary.stat

