#!/bin/bash

inDIR=/mgi_storage/sk/stomics
SN=EXAMPLE-SN
FC=EXAMPLE-FC
outDIR=EXAMPLE-OUT
visualSif=/mgi_storage/sk/stomics/SAW_v5.1.3.sif

export SINGULARITY_BIND=$inDIR,$outDIR

singularity exec ${visualSif} count \
    -i ${outDIR}/00.mapping/${FC}_L01.Aligned.sortedByCoord.out.bam,${outDIR}/00.mapping/${FC}_L02.Aligned.sortedByCoord.out.bam,${outDIR}/00.mapping/${FC}_L03.Aligned.sortedByCoord.out.bam,${outDIR}/00.mapping/${FC}_L04.Aligned.sortedByCoord.out.bam\
    -o ${outDIR}/02.count/${SN}.Aligned.sortedByCoord.out.merge.q10.dedup.target.bam \
    -a ${inDIR}/reference/GRCm39_GENCODE/gencode.vM30.primary_assembly.annotation.gtf \
    -s ${outDIR}/02.count/${SN}.Aligned.sortedByCoord.out.merge.q10.dedup.target.bam.summary.stat \
    -e ${outDIR}/02.count/${SN}.raw.gef \
    --sat_file ${outDIR}/02.count/raw_barcode_gene_exp.txt \
    --umi_on \
    --save_lq \
    --save_dup \
    --sn ${SN} \
    -c 36 \
    -m 128

