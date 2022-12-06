#!/bin/bash

if [[ $# -lt 4 ]] ; then
    echo 'usage: ./_SAW-preWork_G400.sh [SN ID] [FlowCell ID] [output dir] [human/mouse/axolotl/wheat/hamster]'
    echo 'eg. /mgi_storage/sk/luzhigang1/scripts/saw_v5.1_singularity/_SAW-preWork_G400.sh SS200000745BL_C6 V350086216 /mgi_storage/sk/luzhigang1/analysis human'
    exit 1
fi

SN=$1
FC=$2
OUTDIR="$(echo $3 | sed 's/\//\\\//g')"
SPEC=$4

# create SAW output directories
mkdir 00.mapping 01.merge 02.count 03.register 04.tissuecut 05.spatialcluster 06.saturation 07.report analysis 

# mapping script for each lane
if [[ $SPEC = "human" ]]; then
    sed "s/EXAMPLE-SN/${SN}/; s/EXAMPLE-FC/${FC}/; s/EXAMPLE-OUT/${OUTDIR}/; s/GRCm39_GENCODE/GRCh38.p13_GENCODE/" /mgi_storage/sk/luzhigang1/scripts/saw_v5.1_singularity/0_mapping.sh > 00.mapping/s0_mapping_1.sh
    sed "s/EXAMPLE-SN/${SN}/; s/EXAMPLE-FC/${FC}/; s/EXAMPLE-OUT/${OUTDIR}/; s/GRCm39_GENCODE/GRCh38.p13_GENCODE/; s/L01/L02/" /mgi_storage/sk/luzhigang1/scripts/saw_v5.1_singularity/0_mapping.sh > 00.mapping/s0_mapping_2.sh
    sed "s/EXAMPLE-SN/${SN}/; s/EXAMPLE-FC/${FC}/; s/EXAMPLE-OUT/${OUTDIR}/; s/GRCm39_GENCODE/GRCh38.p13_GENCODE/; s/L01/L03/" /mgi_storage/sk/luzhigang1/scripts/saw_v5.1_singularity/0_mapping.sh > 00.mapping/s0_mapping_3.sh
    sed "s/EXAMPLE-SN/${SN}/; s/EXAMPLE-FC/${FC}/; s/EXAMPLE-OUT/${OUTDIR}/; s/GRCm39_GENCODE/GRCh38.p13_GENCODE/; s/L01/L04/" /mgi_storage/sk/luzhigang1/scripts/saw_v5.1_singularity/0_mapping.sh > 00.mapping/s0_mapping_4.sh
elif [[ $SPEC = "wheat" ]]; then
    sed "s/EXAMPLE-SN/${SN}/; s/EXAMPLE-FC/${FC}/; s/EXAMPLE-OUT/${OUTDIR}/; s/GRCm39_GENCODE/IWGSC_RefSeq_v2.1/" /mgi_storage/sk/luzhigang1/scripts/saw_v5.1_singularity/0_mapping.sh > 00.mapping/s0_mapping_1.sh
    sed "s/EXAMPLE-SN/${SN}/; s/EXAMPLE-FC/${FC}/; s/EXAMPLE-OUT/${OUTDIR}/; s/GRCm39_GENCODE/IWGSC_RefSeq_v2.1/; s/L01/L02/" /mgi_storage/sk/luzhigang1/scripts/saw_v5.1_singularity/0_mapping.sh > 00.mapping/s0_mapping_2.sh
    sed "s/EXAMPLE-SN/${SN}/; s/EXAMPLE-FC/${FC}/; s/EXAMPLE-OUT/${OUTDIR}/; s/GRCm39_GENCODE/IWGSC_RefSeq_v2.1/; s/L01/L03/" /mgi_storage/sk/luzhigang1/scripts/saw_v5.1_singularity/0_mapping.sh > 00.mapping/s0_mapping_3.sh
    sed "s/EXAMPLE-SN/${SN}/; s/EXAMPLE-FC/${FC}/; s/EXAMPLE-OUT/${OUTDIR}/; s/GRCm39_GENCODE/IWGSC_RefSeq_v2.1/; s/L01/L04/" /mgi_storage/sk/luzhigang1/scripts/saw_v5.1_singularity/0_mapping.sh > 00.mapping/s0_mapping_4.sh
elif [[ $SPEC = "axolotl" ]]; then
    sed "s/EXAMPLE-SN/${SN}/; s/EXAMPLE-FC/${FC}/; s/EXAMPLE-OUT/${OUTDIR}/; s/GRCm39_GENCODE/axolotl_v6.0_sz/" /mgi_storage/sk/luzhigang1/scripts/saw_v5.1_singularity/0_mapping.sh > 00.mapping/s0_mapping_1.sh
    sed "s/EXAMPLE-SN/${SN}/; s/EXAMPLE-FC/${FC}/; s/EXAMPLE-OUT/${OUTDIR}/; s/GRCm39_GENCODE/axolotl_v6.0_sz/; s/L01/L02/" /mgi_storage/sk/luzhigang1/scripts/saw_v5.1_singularity/0_mapping.sh > 00.mapping/s0_mapping_2.sh
    sed "s/EXAMPLE-SN/${SN}/; s/EXAMPLE-FC/${FC}/; s/EXAMPLE-OUT/${OUTDIR}/; s/GRCm39_GENCODE/axolotl_v6.0_sz/; s/L01/L03/" /mgi_storage/sk/luzhigang1/scripts/saw_v5.1_singularity/0_mapping.sh > 00.mapping/s0_mapping_3.sh
    sed "s/EXAMPLE-SN/${SN}/; s/EXAMPLE-FC/${FC}/; s/EXAMPLE-OUT/${OUTDIR}/; s/GRCm39_GENCODE/axolotl_v6.0_sz/; s/L01/L04/" /mgi_storage/sk/luzhigang1/scripts/saw_v5.1_singularity/0_mapping.sh > 00.mapping/s0_mapping_4.sh
elif [[ $SPEC = "hamster" ]]; then
    sed "s/EXAMPLE-SN/${SN}/; s/EXAMPLE-FC/${FC}/; s/EXAMPLE-OUT/${OUTDIR}/; s/GRCm39_GENCODE/MesAur1.0/" /mgi_storage/sk/luzhigang1/scripts/saw_v5.1_singularity/0_mapping.sh > 00.mapping/s0_mapping_1.sh
    sed "s/EXAMPLE-SN/${SN}/; s/EXAMPLE-FC/${FC}/; s/EXAMPLE-OUT/${OUTDIR}/; s/GRCm39_GENCODE/MesAur1.0/; s/L01/L02/" /mgi_storage/sk/luzhigang1/scripts/saw_v5.1_singularity/0_mapping.sh > 00.mapping/s0_mapping_2.sh
    sed "s/EXAMPLE-SN/${SN}/; s/EXAMPLE-FC/${FC}/; s/EXAMPLE-OUT/${OUTDIR}/; s/GRCm39_GENCODE/MesAur1.0/; s/L01/L03/" /mgi_storage/sk/luzhigang1/scripts/saw_v5.1_singularity/0_mapping.sh > 00.mapping/s0_mapping_3.sh
    sed "s/EXAMPLE-SN/${SN}/; s/EXAMPLE-FC/${FC}/; s/EXAMPLE-OUT/${OUTDIR}/; s/GRCm39_GENCODE/MesAur1.0/; s/L01/L04/" /mgi_storage/sk/luzhigang1/scripts/saw_v5.1_singularity/0_mapping.sh > 00.mapping/s0_mapping_4.sh
else
    sed "s/EXAMPLE-SN/${SN}/; s/EXAMPLE-FC/${FC}/; s/EXAMPLE-OUT/${OUTDIR}/" /mgi_storage/sk/luzhigang1/scripts/saw_v5.1_singularity/0_mapping.sh > 00.mapping/s0_mapping_1.sh
    sed "s/EXAMPLE-SN/${SN}/; s/EXAMPLE-FC/${FC}/; s/EXAMPLE-OUT/${OUTDIR}/; s/L01/L02/" /mgi_storage/sk/luzhigang1/scripts/saw_v5.1_singularity/0_mapping.sh > 00.mapping/s0_mapping_2.sh
    sed "s/EXAMPLE-SN/${SN}/; s/EXAMPLE-FC/${FC}/; s/EXAMPLE-OUT/${OUTDIR}/; s/L01/L03/" /mgi_storage/sk/luzhigang1/scripts/saw_v5.1_singularity/0_mapping.sh > 00.mapping/s0_mapping_3.sh
    sed "s/EXAMPLE-SN/${SN}/; s/EXAMPLE-FC/${FC}/; s/EXAMPLE-OUT/${OUTDIR}/; s/L01/L04/" /mgi_storage/sk/luzhigang1/scripts/saw_v5.1_singularity/0_mapping.sh > 00.mapping/s0_mapping_4.sh
fi

# bcPara files
sed "s/EXAMPLE-SN/${SN}/g; s/EXAMPLE-FC/${FC}/g; s/EXAMPLE-OUT/${OUTDIR}/g" /mgi_storage/sk/luzhigang1/scripts/saw_v5.1_singularity/G400.bcPara > 00.mapping/L01.bcPara 
sed "s/EXAMPLE-SN/${SN}/g; s/EXAMPLE-FC/${FC}/g; s/EXAMPLE-OUT/${OUTDIR}/g; s/L01/L02/g" /mgi_storage/sk/luzhigang1/scripts/saw_v5.1_singularity/G400.bcPara > 00.mapping/L02.bcPara 
sed "s/EXAMPLE-SN/${SN}/g; s/EXAMPLE-FC/${FC}/g; s/EXAMPLE-OUT/${OUTDIR}/g; s/L01/L03/g" /mgi_storage/sk/luzhigang1/scripts/saw_v5.1_singularity/G400.bcPara > 00.mapping/L03.bcPara 
sed "s/EXAMPLE-SN/${SN}/g; s/EXAMPLE-FC/${FC}/g; s/EXAMPLE-OUT/${OUTDIR}/g; s/L01/L04/g" /mgi_storage/sk/luzhigang1/scripts/saw_v5.1_singularity/G400.bcPara > 00.mapping/L04.bcPara 

# merge script
sed "s/EXAMPLE-SN/${SN}/; s/EXAMPLE-FC/${FC}/; s/EXAMPLE-OUT/${OUTDIR}/" /mgi_storage/sk/luzhigang1/scripts/saw_v5.1_singularity/1_merge.sh > 01.merge/s1_merge.sh

# count script
if [[ $SPEC = "human" ]]; then
    sed "s/EXAMPLE-SN/${SN}/; s/EXAMPLE-FC/${FC}/; s/EXAMPLE-OUT/${OUTDIR}/; s/GRCm39_GENCODE/GRCh38.p13_GENCODE/; s/vM30/v41/" /mgi_storage/sk/luzhigang1/scripts/saw_v5.1_singularity/2_count.sh > 02.count/s2_count.sh
elif [[ $SPEC = "wheat" ]]; then
    sed "s/EXAMPLE-SN/${SN}/; s/EXAMPLE-FC/${FC}/; s/EXAMPLE-OUT/${OUTDIR}/; s/GRCm39_GENCODE/IWGSC_RefSeq_v2.1/; s/gencode.vM30.primary_assembly.annotation.gtf/iwgsc_refseqv2.1_HC-LC_chloro_mito_corr.gff3/" /mgi_storage/sk/luzhigang1/scripts/saw_v5.1_singularity/2_count.sh > 02.count/s2_count.sh
elif [[ $SPEC = "axolotl" ]]; then
    sed "s/EXAMPLE-SN/${SN}/; s/EXAMPLE-FC/${FC}/; s/EXAMPLE-OUT/${OUTDIR}/; s/GRCm39_GENCODE/axolotl_v6.0_sz/; s/gencode.vM30.primary_assembly.annotation.gtf/AmexG_v6.0_sz_ed_corr.gtf/" /mgi_storage/sk/luzhigang1/scripts/saw_v5.1_singularity/2_count.sh > 02.count/s2_count.sh
elif [[ $SPEC = "hamster" ]]; then
    sed "s/EXAMPLE-SN/${SN}/; s/EXAMPLE-FC/${FC}/; s/EXAMPLE-OUT/${OUTDIR}/; s/GRCm39_GENCODE/MesAur1.0/; s/gencode.vM30.primary_assembly.annotation.gtf/Mesocricetus_auratus.MesAur1.0.107.gtf/" /mgi_storage/sk/luzhigang1/scripts/saw_v5.1_singularity/2_count.sh > 02.count/s2_count.sh
else
    sed "s/EXAMPLE-SN/${SN}/; s/EXAMPLE-FC/${FC}/; s/EXAMPLE-OUT/${OUTDIR}/" /mgi_storage/sk/luzhigang1/scripts/saw_v5.1_singularity/2_count.sh > 02.count/s2_count.sh
fi

# register
sed "s/EXAMPLE-SN/${SN}/; s/EXAMPLE-OUT/${OUTDIR}/" /mgi_storage/sk/luzhigang1/scripts/saw_v5.1_singularity/3_register.sh > 03.register/s3_register.sh
sed "s/EXAMPLE-SN/${SN}/; s/EXAMPLE-OUT/${OUTDIR}/" /mgi_storage/sk/luzhigang1/scripts/saw_v5.1_singularity/3_rapidRegister.sh > 03.register/s3_rapidRegister.sh
sed "s/EXAMPLE-SN/${SN}/; s/EXAMPLE-OUT/${OUTDIR}/" /mgi_storage/sk/luzhigang1/scripts/saw_v5.1_singularity/3_ipr2img.sh > 03.register/s3_ipr2img.sh

# tissueCut and cellCut script
sed "s/EXAMPLE-SN/${SN}/; s/EXAMPLE-OUT/${OUTDIR}/" /mgi_storage/sk/luzhigang1/scripts/saw_v5.1_singularity/4_tissueCut.sh > 04.tissuecut/s4_tissueCut.sh
sed "s/EXAMPLE-SN/${SN}/; s/EXAMPLE-OUT/${OUTDIR}/" /mgi_storage/sk/luzhigang1/scripts/saw_v5.1_singularity/4_cellCut.sh > 04.tissuecut/s4_cellCut.sh

# spatialCluster and cellCluster script
sed "s/EXAMPLE-SN/${SN}/; s/EXAMPLE-OUT/${OUTDIR}/" /mgi_storage/sk/luzhigang1/scripts/saw_v5.1_singularity/5_spatialCluster.sh > 05.spatialcluster/s5_spatialCluster.sh
sed "s/EXAMPLE-SN/${SN}/; s/EXAMPLE-OUT/${OUTDIR}/" /mgi_storage/sk/luzhigang1/scripts/saw_v5.1_singularity/5_cellCluster.sh > 05.spatialcluster/s5_cellCluster.sh

# saturation script
sed "s/EXAMPLE-SN/${SN}/; s/EXAMPLE-FC/${FC}/; s/EXAMPLE-OUT/${OUTDIR}/" /mgi_storage/sk/luzhigang1/scripts/saw_v5.1_singularity/6_saturation.sh > 06.saturation/s6_saturation.sh

# report script
sed "s/EXAMPLE-SN/${SN}/; s/EXAMPLE-FC/${FC}/; s/EXAMPLE-OUT/${OUTDIR}/; s/EXAMPLE-SPEC/${SPEC}/" /mgi_storage/sk/luzhigang1/scripts/saw_v5.1_singularity/7_report.sh > 07.report/s7_report.sh

# other scripts
sed "s/EXAMPLE-SN/${SN}/; s/EXAMPLE-OUT/${OUTDIR}/" /mgi_storage/sk/luzhigang1/scripts/saw_v5.1_singularity/CIDCount.sh > CIDCount.sh
sed "s/EXAMPLE-SN/${SN}/; s/EXAMPLE-OUT/${OUTDIR}/" /mgi_storage/sk/luzhigang1/scripts/saw_v5.1_singularity/cellCut_tools.sh > cellCut_tools.sh

# Seurat script
#sed "s/EXAMPLE-SN/${SN}/; s/EXAMPLE-OUT/${OUTDIR}/; s/EXAMPLE-SPEC/${SPEC}/" /mgi_storage/sk/luzhigang1/scripts/Bin_QC_Seurat.sh > analysis/Bin_QC_Seurat.sh
sed "s/EXAMPLE-SN/${SN}/; s/EXAMPLE-OUT/${OUTDIR}/; s/EXAMPLE-SPEC/${SPEC}/" /mgi_storage/sk/luzhigang1/scripts/run_SCT_BayesSpace_Giotto.sh > analysis/run_SCT_BayesSpace_Giotto.sh

echo "Scripts copied. Here are the input dirs/files:"
find /mgi_storage/sk/stomics/ -type d -name ${FC}
find /mgi_storage/sk/stomics/ -type f -name ${SN}*.h5
find /mgi_storage/sk/stomics/ -type d -name ${SN}

echo " "
echo "Current pipeline: mapping->merge->count->rapidRegister+ipr2img->tissueCut->spatialCluster->saturation->report"
