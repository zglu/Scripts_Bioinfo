#!/bin/bash

if [[ $# -lt 5 ]] ; then
    echo 'usage: ./_SAW-preWork_T7.sh [SN ID] [FlowCell ID] [output dir] [fastq dir] [sample barcode] [human/mouse]'
    echo 'eg. ./_SAW-preWork_T7.sh SS200000XXXTR_D4 E100026XX0 /mgi_storage/sk/luzhigang1/analysis /mgi_storage/sk/stomics/fastq/E100026460 12 mouse'
    exit 1
fi

SN=$1
FC=$2
OUTDIR="$(echo $3 | sed 's/\//\\\//g')"

if [[ $4 = "0" ]]; then
    FQDIR="$(echo "/mgi_storage/sk/stomics/fastq/${FC}" | sed 's/\//\\\//g')"
else
    FQDIR="$(echo $4 | sed 's/\//\\\//g')"
fi

SB=$5
SPEC=$6

# create SAW output directories
mkdir 00.mapping 01.merge 02.count 03.register 04.tissuecut 05.spatialcluster 06.saturation 07.report

# mapping script
if [[ $6 = "human" ]]; then
    sed "s/EXAMPLE-SN/${SN}/; s/EXAMPLE-FC/${FC}/; s/EXAMPLE-OUT/${OUTDIR}/; s/EXAMPLE-FQ/${FQDIR}/; s/GRCm39_GENCODE/GRCh38.p13_GENCODE/" /mgi_storage/sk/luzhigang1/scripts/saw/0_mapping_T7.sh > 00.mapping/s0_mapping_T7.sh
else
    sed "s/EXAMPLE-SN/${SN}/; s/EXAMPLE-FC/${FC}/; s/EXAMPLE-OUT/${OUTDIR}/; s/EXAMPLE-FQ/${FQDIR}/" /mgi_storage/sk/luzhigang1/scripts/saw/0_mapping_T7.sh > 00.mapping/s0_mapping_T7.sh
fi

# bcPara files
sed "s/EXAMPLE-SN/${SN}/; s/EXAMPLE-FC/${FC}/; s/EXAMPLE-OUT/${OUTDIR}/; s/EXAMPLE-FQ/${FQDIR}/; s/EXAMPLE-SB/${SB}/" /mgi_storage/sk/luzhigang1/scripts/saw/T7.bcPara > 00.mapping/L01.bcPara 

# count script
if [[ $6 = "human" ]]; then
    sed "s/EXAMPLE-SN/${SN}/; s/EXAMPLE-FC/${FC}/; s/EXAMPLE-OUT/${OUTDIR}/; s/GRCm39_GENCODE/GRCh38.p13_GENCODE/; s/vM30/v41/" /mgi_storage/sk/luzhigang1/scripts/saw/2_count_T7.sh > 02.count/s2_count_T7.sh
else
    sed "s/EXAMPLE-SN/${SN}/; s/EXAMPLE-FC/${FC}/; s/EXAMPLE-OUT/${OUTDIR}/" /mgi_storage/sk/luzhigang1/scripts/saw/2_count_T7.sh > 02.count/s2_count_T7.sh
fi

# tissueCut script
sed "s/EXAMPLE-SN/${SN}/; s/EXAMPLE-OUT/${OUTDIR}/" /mgi_storage/sk/luzhigang1/scripts/saw/4_tissueCut_T7.sh > 04.tissuecut/s4_tissueCut_T7.sh

# spatialCluster script
sed "s/EXAMPLE-SN/${SN}/; s/EXAMPLE-OUT/${OUTDIR}/" /mgi_storage/sk/luzhigang1/scripts/saw/5_spatialCluster.sh > 05.spatialcluster/s5_spatialCluster.sh

# saturation script
sed "s/EXAMPLE-SN/${SN}/; s/EXAMPLE-FC/${FC}/; s/EXAMPLE-OUT/${OUTDIR}/" /mgi_storage/sk/luzhigang1/scripts/saw/6_saturation_T7.sh > 06.saturation/s6_saturation_T7.sh

# report script
sed "s/EXAMPLE-SN/${SN}/; s/EXAMPLE-FC/${FC}/; s/EXAMPLE-OUT/${OUTDIR}/" /mgi_storage/sk/luzhigang1/scripts/saw/7_report_T7.sh > 07.report/s7_report_T7.sh

# other scripts
sed "s/EXAMPLE-SN/${SN}/; s/EXAMPLE-OUT/${OUTDIR}/" /mgi_storage/sk/luzhigang1/scripts/saw/CIDcount.sh > CIDcount.sh
sed "s/EXAMPLE-SN/${SN}/; s/EXAMPLE-OUT/${OUTDIR}/" /mgi_storage/sk/luzhigang1/scripts/saw/gefTools.sh > gefTools.sh

echo "Scripts copied."
echo "Please check reference genome in s0_mapping.sh and s2_count.sh"
