#!/bin/bash

inDIR=/mgi_storage/sk/luzhigang1/scripts
outDIR=EXAMPLE-OUT
SN=EXAMPLE-SN
SPEC=EXAMPLE-SPEC
SeuratSif=/mgi_storage/sk/stomics/Seurat_4.1.0.sif

export SINGULARITY_BIND=$inDIR,$outDIR

singularity exec ${SeuratSif} Rscript \
    ${inDIR}/a1_bin2rds_distr.R ${outDIR}/Bin_Analysis/${SN}.raw.gem 50 ${SPEC}_bin50

singularity exec ${SeuratSif} Rscript \
    ${inDIR}/a1_bin2rds_distr.R ${outDIR}/Bin_Analysis/${SN}.raw.gem 200 ${SPEC}_bin200

singularity exec ${SeuratSif} Rscript \
    ${inDIR}/a0_bin2rds-SCT.R ${outDIR}/Bin_Analysis/${SN}.tissue.bin1.gem 50 ${SPEC}_bin50

#cellbin
#singularity exec ${SeuratSif} Rscript \
#    ${inDIR}/b0_cellbin2rds_sct.R ${outDIR}/Bin_Analysis/${SN}.cellbin.gem ${SPEC}

# check feaures if provided
#for i in `cat markers.ids`; do
#  echo $i
#  singularity exec ${SeuratSif} Rscript ${inDIR}/featureSpatial_raw.R ${SN}.raw.gem_BIN200.rds $i
#  singularity exec ${SeuratSif} Rscript ${inDIR}/featureSpatial_sct.R ${SN}.tissue.bin1.gem_BIN50.rds_SCT.rds $i
#done

