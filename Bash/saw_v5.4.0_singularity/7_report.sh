#!/bin/bash

inDIR=/mgi_storage/sk/stomics
outDIR=EXAMPLE-OUT
SN=EXAMPLE-SN
FC=EXAMPLE-FC
visualSif=/mgi_storage/sk/stomics/SAW_v5.4.0.sif

imageIPR=$(find ${outDIR}/03.register -maxdepth 1 -name ${SN}*.ipr | head -1)

export SINGULARITY_BIND=$inDIR,$outDIR

singularity exec ${visualSif} cellCut bgef \
    -i ${outDIR}/02.count/${SN}.raw.gef \
    -o ${outDIR}/04.tissuecut/${SN}.gef \
    -b 1,10,20,50,100,200,500 \
    -O Transcriptomics

singularity exec ${visualSif} report \
    -m ${outDIR}/00.mapping/${FC}_L01_barcodeMap.stat,${outDIR}/00.mapping/${FC}_L02_barcodeMap.stat,${outDIR}/00.mapping/${FC}_L03_barcodeMap.stat,${outDIR}/00.mapping/${FC}_L04_barcodeMap.stat \
    -a ${outDIR}/00.mapping/${FC}_L01.Log.final.out,${outDIR}/00.mapping/${FC}_L02.Log.final.out,${outDIR}/00.mapping/${FC}_L03.Log.final.out,${outDIR}/00.mapping/${FC}_L04.Log.final.out \
    -g ${outDIR}/02.count/${SN}.Aligned.sortedByCoord.out.merge.q10.dedup.target.bam.summary.stat \
    -l ${outDIR}/04.tissuecut/tissuecut.stat \
    -n ${outDIR}/04.tissuecut/${SN}.gef \
    -d ${outDIR}/05.spatialcluster/${SN}.spatial.cluster.h5ad \
    -t ${outDIR}/06.saturation/plot_200x200_saturation.png \
    -b ${outDIR}/04.tissuecut/tissue_fig/scatter_200x200_MID_gene_counts.png \
    -v ${outDIR}/04.tissuecut/tissue_fig/violin_200x200_MID_gene.png \
    -c ${outDIR}/04.tissuecut/tissue_fig/statistic_200x200_MID_gene_DNB.png \
    --bin50Saturation ${outDIR}/04.tissuecut/tissue_fig/scatter_50x50_MID_gene_counts.png \
    --bin50violin ${outDIR}/04.tissuecut/tissue_fig/violin_50x50_MID_gene.png \
    --bin50MIDGeneDNB ${outDIR}/04.tissuecut/tissue_fig/statistic_50x50_MID_gene_DNB.png \
    --bin100Saturation ${outDIR}/04.tissuecut/tissue_fig/scatter_100x100_MID_gene_counts.png \
    --bin100violin ${outDIR}/04.tissuecut/tissue_fig/violin_100x100_MID_gene.png \
    --bin100MIDGeneDNB ${outDIR}/04.tissuecut/tissue_fig/statistic_100x100_MID_gene_DNB.png \
    --bin150Saturation ${outDIR}/04.tissuecut/tissue_fig/scatter_150x150_MID_gene_counts.png \
    --bin150violin ${outDIR}/04.tissuecut/tissue_fig/violin_150x150_MID_gene.png \
    --bin150MIDGeneDNB ${outDIR}/04.tissuecut/tissue_fig/statistic_150x150_MID_gene_DNB.png \
    -r standard_version \
    --pipelineVersion SAW_v5.4.0 \
    -i ${outDIR}/03.register/${SN}.rpi \
    --iprFile ${imageIPR} \
    -s ${SN} \
    -o ${outDIR}/07.report \
    --species EXAMPLE-SPEC\
    --tissue tissue \
    --reference ref_source

# delete lines imageIPR... -i... --iprFile... if no register/rapidRegister was done

# add the following lines if register/cellcut/cellcluster were done
#    --cellBinGef ${outDIR}/04.cellcut/${SN}.cellbin.gef \
#    --cellCluster ${outDIR}/05.cellcluster/${SN}.cell.cluster.h5ad \