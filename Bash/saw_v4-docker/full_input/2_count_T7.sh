#!/bin/bash

SN=EXAMPLE-SN
FC=EXAMPLE-FC
outDIR=EXAMPLE-OUT
refDIR=/mgi_storage/sk/stomics/reference/GRCm39_GENCODE

docker run --rm --cpus=64 --memory=130g -u $(id -u):$(id -g) \
  -v ${outDIR}:/outDir \
  -v ${refDIR}:/refDir \
stomics/saw:04.1.0 /bin/bash count \
    -i /outDir/00.mapping/${FC}.Aligned.sortedByCoord.out.bam\
    -o /outDir/02.count/${SN}.Aligned.sortedByCoord.out.merge.q10.dedup.target.bam \
    -a /refDir/gencode.vM30.primary_assembly.annotation.gtf \
    -s /outDir/02.count/${SN}.Aligned.sortedByCoord.out.merge.q10.dedup.target.bam.summary.stat \
    -e /outDir/02.count/${SN}.raw.gef \
    --sat_file /outDir/02.count/raw_barcode_gene_exp.txt \
    --umi_on \
    --save_lq \
    --save_dup \
    --sn ${SN} \
    -c 36 \
    -m 128

