#!/bin/bash

inDIR=/mgi_storage/sk/stomics
SN=SS200000183TR_D4
FC1=E100043446
FC2=E100042432
FC3=V300109183
FC4=V300109651
outDIR=EXAMPLE-OUT
visualSif=/mgi_storage/sk/stomics/SAW_v5.1.3.sif

export SINGULARITY_BIND=$inDIR,$outDIR

singularity exec ${visualSif} count \
    -i ${outDIR}/${FC1}/00.mapping/${FC1}.Aligned.sortedByCoord.out.bam,${outDIR}/${FC2}/00.mapping/${FC2}.Aligned.sortedByCoord.out.bam,${outDIR}/${FC3}/00.mapping/${FC3}_L01.Aligned.sortedByCoord.out.bam,${outDIR}/${FC3}/00.mapping/${FC3}_L02.Aligned.sortedByCoord.out.bam,${outDIR}/${FC3}/00.mapping/${FC3}_L03.Aligned.sortedByCoord.out.bam,${outDIR}/${FC3}/00.mapping/${FC3}_L04.Aligned.sortedByCoord.out.bam,${outDIR}/${FC4}/00.mapping/${FC4}_L01.Aligned.sortedByCoord.out.bam,${outDIR}/${FC4}/00.mapping/${FC4}_L02.Aligned.sortedByCoord.out.bam,${outDIR}/${FC4}/00.mapping/${FC4}_L03.Aligned.sortedByCoord.out.bam,${outDIR}/${FC4}/00.mapping/${FC4}_L04.Aligned.sortedByCoord.out.bam \
    -o ${outDIR}/02.count/${SN}.Aligned.sortedByCoord.out.merge.q10.dedup.target.bam \
    -a ${inDIR}/reference/GRCm39_GENCODE/gencode.vM30.primary_assembly.annotation.gtf \
    -s ${outDIR}/02.count/${SN}.Aligned.sortedByCoord.out.merge.q10.dedup.target.bam.summary.stat \
    -e ${outDIR}/02.count/${SN}.raw.gef \
    --sat_file ${outDIR}/02.count/raw_barcode_gene_exp.txt \
    --umi_on \
    --save_lq \
    --save_dup \
    --sn ${SN} \
    -c 64 \
    -m 320

