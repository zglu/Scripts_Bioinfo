#!/bin/bash

if [[ $# -lt 4 ]] ; then
    echo 'usage: ./_SAW-preWork_T7.sh [SN ID] [FlowCell ID] [output dir] [human/mouse]'
    echo 'eg. ./_SAW-preWork_T7.sh SS200000XXXTR_D4 E100026XX0 /mgi_storage/sk/luzhigang1/analysis mouse'
    exit 1
fi

SN=$1
FC=$2
OUTDIR="$(echo $3 | sed 's/\//\\\//g')"
SPEC=$4

# create SAW output directories
mkdir 00.mapping 01.merge 02.count 03.register 04.tissuecut 05.spatialcluster 06.saturation 07.report

# mapping script
if [[ $SPEC = "human" ]]; then
    sed "s/EXAMPLE-SN/${SN}/; s/EXAMPLE-FC/${FC}/; s/EXAMPLE-OUT/${OUTDIR}/; s/GRCm39_GENCODE/GRCh38.p13_GENCODE/" ../0_mapping_T7.sh > 00.mapping/s0_mapping_T7.sh
elif [[ $SPEC = "wheet" ]]; then
    sed "s/EXAMPLE-SN/${SN}/; s/EXAMPLE-FC/${FC}/; s/EXAMPLE-OUT/${OUTDIR}/; s/GRCm39_GENCODE/Triticum_aestivum.IWGSC54/" ../0_mapping_T7.sh > 00.mapping/s0_mapping_T7.sh
elif [[ $SPEC = "axolotl" ]]; then
    sed "s/EXAMPLE-SN/${SN}/; s/EXAMPLE-FC/${FC}/; s/EXAMPLE-OUT/${OUTDIR}/; s/GRCm39_GENCODE/axolotl_v6.0_sz/" ../0_mapping_T7.sh > 00.mapping/s0_mapping_T7.sh
elif [[ $SPEC = "hamster" ]]; then
    sed "s/EXAMPLE-SN/${SN}/; s/EXAMPLE-FC/${FC}/; s/EXAMPLE-OUT/${OUTDIR}/; s/GRCm39_GENCODE/MesAur1.0/" ../0_mapping_T7.sh > 00.mapping/s0_mapping_T7.sh
else
    sed "s/EXAMPLE-SN/${SN}/; s/EXAMPLE-FC/${FC}/; s/EXAMPLE-OUT/${OUTDIR}/" ../0_mapping_T7.sh > 00.mapping/s0_mapping_T7.sh 
fi

# bcPara files
sed "s/EXAMPLE-SN/${SN}/g; s/EXAMPLE-FC/${FC}/g; s/EXAMPLE-OUT/${OUTDIR}/g" ../T7.bcPara > 00.mapping/L01.bcPara 

# count script
if [[ $SPEC = "human" ]]; then
    sed "s/EXAMPLE-SN/${SN}/; s/EXAMPLE-FC/${FC}/; s/EXAMPLE-OUT/${OUTDIR}/; s/GRCm39_GENCODE/GRCh38.p13_GENCODE/; s/vM30/v41/" ../2_count_T7.sh > 02.count/s2_count_T7.sh
elif [[ $SPEC = "wheet" ]]; then
    sed "s/EXAMPLE-SN/${SN}/; s/EXAMPLE-FC/${FC}/; s/EXAMPLE-OUT/${OUTDIR}/; s/GRCm39_GENCODE/Triticum_aestivum.IWGSC54/; s/gencode.vM30.primary_assembly.annotation.gtf/Triticum_aestivum.IWGSC.54.gff3/" ../2_count_T7.sh > 02.count/s2_count_T7.sh
elif [[ $SPEC = "axolotl" ]]; then
    sed "s/EXAMPLE-SN/${SN}/; s/EXAMPLE-FC/${FC}/; s/EXAMPLE-OUT/${OUTDIR}/; s/GRCm39_GENCODE/axolotl_v6.0_sz/; s/gencode.vM30.primary_assembly.annotation.gtf/AmexT_v47.gtf/" ../2_count_T7.sh > 02.count/s2_count_T7.sh
elif [[ $SPEC = "hamster" ]]; then
    sed "s/EXAMPLE-SN/${SN}/; s/EXAMPLE-FC/${FC}/; s/EXAMPLE-OUT/${OUTDIR}/; s/GRCm39_GENCODE/MesAur1.0/; s/gencode.vM30.primary_assembly.annotation.gtf/Mesocricetus_auratus.MesAur1.0.107.gtf/" ../2_count_T7.sh > 02.count/s2_count_T7.sh
else
    sed "s/EXAMPLE-SN/${SN}/; s/EXAMPLE-FC/${FC}/; s/EXAMPLE-OUT/${OUTDIR}/" ../2_count_T7.sh > 02.count/s2_count_T7.sh
fi

# register
sed "s/EXAMPLE-SN/${SN}/; s/EXAMPLE-OUT/${OUTDIR}/" ../3_register.sh > 03.register/s3_register.sh
sed "s/EXAMPLE-SN/${SN}/; s/EXAMPLE-OUT/${OUTDIR}/" ../3_rapidRegister.sh > 03.register/s3_rapidRegister.sh
sed "s/EXAMPLE-SN/${SN}/; s/EXAMPLE-OUT/${OUTDIR}/" ../3_ipr2img.sh > 03.register/s3_ipr2img.sh

# tissueCut script
sed "s/EXAMPLE-SN/${SN}/; s/EXAMPLE-OUT/${OUTDIR}/; s/EXAMPLE-FC/${FC}/" ../4_tissueCut_T7.sh > 04.tissuecut/s4_tissueCut_T7.sh
sed "s/EXAMPLE-SN/${SN}/; s/EXAMPLE-OUT/${OUTDIR}/" ../4_cellCut.sh > 04.tissuecut/s4_cellCut.sh

# spatialCluster script
sed "s/EXAMPLE-SN/${SN}/; s/EXAMPLE-OUT/${OUTDIR}/" ../5_spatialCluster.sh > 05.spatialcluster/s5_spatialCluster.sh
sed "s/EXAMPLE-SN/${SN}/; s/EXAMPLE-OUT/${OUTDIR}/" ../5_cellCluster.sh > 05.spatialcluster/s5_cellCluster.sh

# saturation script
sed "s/EXAMPLE-SN/${SN}/; s/EXAMPLE-FC/${FC}/; s/EXAMPLE-OUT/${OUTDIR}/" ../6_saturation_T7.sh > 06.saturation/s6_saturation_T7.sh

# report script
sed "s/EXAMPLE-SN/${SN}/; s/EXAMPLE-FC/${FC}/; s/EXAMPLE-OUT/${OUTDIR}/; s/EXAMPLE-SPEC/${SPEC}/" ../7_report_T7.sh > 07.report/s7_report_T7.sh

# other scripts
sed "s/EXAMPLE-SN/${SN}/; s/EXAMPLE-OUT/${OUTDIR}/" ../CIDcount.sh > CIDcount.sh
sed "s/EXAMPLE-SN/${SN}/; s/EXAMPLE-OUT/${OUTDIR}/" ../cellCut_tools.sh > cellCut_tools.sh
