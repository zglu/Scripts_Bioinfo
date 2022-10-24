#!/bin/bash

SN=EXAMPLE-SN
FC=EXAMPLE-FC
outDIR=EXAMPLE-OUT

docker run --rm --cpus=12 --memory=20g \
  -v ${outDIR}:/outDir \
stomics/saw:04.1.0 /bin/bash report \
    -m /outDir/00.mapping/${FC}_barcodeMap.stat \
    -a /outDir/00.mapping/${FC}.Log.final.out \
    -g /outDir/02.count/${SN}.Aligned.sortedByCoord.out.merge.q10.dedup.target.bam.summary.stat \
    -l /outDir/04.tissuecut/tissuecut.stat \
    -n /outDir/04.tissuecut/${SN}.gef \
    -d /outDir/05.spatialcluster/${SN}.spatial.cluster.h5ad \
    -t /outDir/06.saturation/plot_200x200_saturation.png \
    -b /outDir/04.tissuecut/tissue_fig/scatter_200x200_MID_gene_counts.png \
    -v /outDir/04.tissuecut/tissue_fig/violin_200x200_MID_gene.png \
#    -i /outDir/04.tissuecut/tissue_fig/${SN}.ssDNA.rpi \
    -o /outDir/07.report \
    -r standard_version \
    --pipelineVersion SAW_v4.1.0 \
    -s ${SN} &&\

