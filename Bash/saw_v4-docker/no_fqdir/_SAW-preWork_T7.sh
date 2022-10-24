#!/bin/bash

if [[ $# -lt 4 ]] ; then
    echo 'usage: ./_SAW-preWork_T7.sh [SN ID] [FlowCell ID] [output dir] [sample barcode] [human/mouse]'
    echo 'eg. ./_SAW-preWork_T7.sh SS200000XXXTR_D4 E100026XX0 /mgi_storage/sk/luzhigang1/analysis 12 mouse'
    exit 1
fi

SN=$1
FC=$2
OUTDIR="$(echo $3 | sed 's/\//\\\//g')"
SB=$4
SPEC=$5

# create SAW output directories
mkdir 00.mapping 01.merge 02.count 03.register 04.tissuecut 05.spatialcluster 06.saturation 07.report

# mapping script
if [[ $5 = "human" ]]; then
    sed "s/EXAMPLE-SN/${SN}/; s/EXAMPLE-FC/${FC}/; s/EXAMPLE-OUT/${OUTDIR}/; s/GRCm39_GENCODE/GRCh38.p13_GENCODE/" /mgi_storage/sk/luzhigang1/scripts/saw/0_mapping_T7.sh > 00.mapping/0_mapping_T7.sh
else
    sed "s/EXAMPLE-SN/${SN}/; s/EXAMPLE-FC/${FC}/; s/EXAMPLE-OUT/${OUTDIR}/" /mgi_storage/sk/luzhigang1/scripts/saw/0_mapping_T7.sh > 00.mapping/0_mapping_T7.sh
fi

# bcPara files
sed "s/EXAMPLE-SN/${SN}/; s/EXAMPLE-FC/${FC}/; s/EXAMPLE-OUT/${OUTDIR}/; s/EXAMPLE-SB/${SB}/" /mgi_storage/sk/luzhigang1/scripts/saw/T7.bcPara > 00.mapping/L01.bcPara 

# count script
if [[ $5 = "human" ]]; then
    sed "s/EXAMPLE-SN/${SN}/; s/EXAMPLE-FC/${FC}/; s/EXAMPLE-OUT/${OUTDIR}/; s/GRCm39_GENCODE/GRCh38.p13_GENCODE/; s/vM30/v41/" /mgi_storage/sk/luzhigang1/scripts/saw/2_count_T7.sh > 02.count/2_count_T7.sh
else
    sed "s/EXAMPLE-SN/${SN}/; s/EXAMPLE-FC/${FC}/; s/EXAMPLE-OUT/${OUTDIR}/" /mgi_storage/sk/luzhigang1/scripts/saw/2_count_T7.sh > 02.count/2_count_T7.sh
fi

# tissueCut script
sed "s/EXAMPLE-SN/${SN}/; s/EXAMPLE-FC/${FC}/; s/EXAMPLE-OUT/${OUTDIR}/" /mgi_storage/sk/luzhigang1/scripts/saw/4_tissueCut_T7.sh > 04.tissuecut/4_tissueCut_T7.sh

# spatialCluster script
sed "s/EXAMPLE-SN/${SN}/; s/EXAMPLE-OUT/${OUTDIR}/" /mgi_storage/sk/luzhigang1/scripts/saw/5_spatialCluster.sh > 05.spatialcluster/5_spatialCluster.sh

# saturation script
sed "s/EXAMPLE-SN/${SN}/; s/EXAMPLE-FC/${FC}/; s/EXAMPLE-OUT/${OUTDIR}/" /mgi_storage/sk/luzhigang1/scripts/saw/6_saturation_T7.sh > 06.saturation/6_saturation_T7.sh

# report script
sed "s/EXAMPLE-SN/${SN}/; s/EXAMPLE-FC/${FC}/; s/EXAMPLE-OUT/${OUTDIR}/" /mgi_storage/sk/luzhigang1/scripts/saw/7_report_T7.sh > 07.report/7_report_T7.sh

# other scripts
sed "s/EXAMPLE-SN/${SN}/; s/EXAMPLE-OUT/${OUTDIR}/" /mgi_storage/sk/luzhigang1/scripts/saw/CIDcount.sh > CIDcount.sh
sed "s/EXAMPLE-SN/${SN}/; s/EXAMPLE-OUT/${OUTDIR}/" /mgi_storage/sk/luzhigang1/scripts/saw/gefTools.sh > gefTools.sh

echo "Please check reference genome in 0_mapping.sh and 2_count.sh"
