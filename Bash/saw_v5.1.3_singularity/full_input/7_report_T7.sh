#!/bin/bash

visualSif=/mgi_storage/sk/stomics/SAW_v5.1.3.sif
SN=EXAMPLE-SN
FC=EXAMPLE-FC
outDIR=EXAMPLE-OUT

export SINGULARITY_BIND=$outDIR

singularity exec ${visualSif} report \
    -m ${outDIR}/00.mapping/${FC}_barcodeMap.stat \
    -a ${outDIR}/00.mapping/${FC}.Log.final.out \
    -g ${outDIR}/02.count/${SN}.Aligned.sortedByCoord.out.merge.q10.dedup.target.bam.summary.stat \
    -l ${outDIR}/04.tissuecut/tissuecut.stat \
    -n ${outDIR}/04.tissuecut/${SN}.gef \
    -d ${outDIR}/05.spatialcluster/${SN}.spatial.cluster.h5ad \
    -t ${outDIR}/06.saturation/plot_200x200_saturation.png \
    -b ${outDIR}/04.tissuecut/tissue_fig/scatter_200x200_MID_gene_counts.png \
    -v ${outDIR}/04.tissuecut/tissue_fig/violin_200x200_MID_gene.png \
#    -i ${outDIR}/04.tissuecut/tissue_fig/${SN}.ssDNA.rpi \
    -o ${outDIR}/07.report \
    -r standard_version \
    --pipelineVersion SAW_v4.1.0 \
    -s ${SN}

