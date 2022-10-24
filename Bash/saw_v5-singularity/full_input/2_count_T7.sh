#!/bin/bash

visualSif=/mgi_storage/sk/stomics/SAW_v5.1.3.sif
SN=EXAMPLE-SN
FC=EXAMPLE-FC
outDIR=EXAMPLE-OUT
refDIR=/mgi_storage/sk/stomics/reference/GRCm39_GENCODE

export SINGULARITY_BIND=$outDIR,$refDIR

singularity exec ${visualSif} count \
    -i ${outDIR}/00.mapping/${FC}.Aligned.sortedByCoord.out.bam\
    -o ${outDIR}/02.count/${SN}.Aligned.sortedByCoord.out.merge.q10.dedup.target.bam \
    -a ${refDIR}/gencode.vM30.primary_assembly.annotation.gtf \
    -s ${outDIR}/02.count/${SN}.Aligned.sortedByCoord.out.merge.q10.dedup.target.bam.summary.stat \
    -e ${outDIR}/02.count/${SN}.raw.gef \
    --sat_file ${outDIR}/02.count/raw_barcode_gene_exp.txt \
    --umi_on \
    --save_lq \
    --save_dup \
    --sn ${SN} \
    -c 36 \
    -m 128

