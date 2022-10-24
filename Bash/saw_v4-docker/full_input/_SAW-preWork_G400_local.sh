#!/bin/bash

if [[ $# -lt 4 ]] ; then
    echo 'usage: ./_SAW-preWork_G400.sh [SN ID] [FlowCell ID] [output dir] [fastq dir] [human/mouse]'
    echo 'eg. /mgi_storage/sk/luzhigang1/scripts/saw/_SAW-preWork_G400.sh SS200000745BL_C6 V350086216 /mgi_storage/sk/luzhigang1/analysis /mgi_storage/sk/stomics/fastq/V350086216 human'
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

SPEC=$5

# create SAW output directories
mkdir 00.mapping 01.merge 02.count 03.register 04.tissuecut 05.spatialcluster 06.saturation 07.report

# mapping script for each lane
if [[ $5 = "human" ]]; then
    sed "s/EXAMPLE-SN/${SN}/; s/EXAMPLE-FC/${FC}/; s/EXAMPLE-OUT/${OUTDIR}/; s/EXAMPLE-FQ/${FQDIR}/; s/GRCm39_GENCODE/GRCh38.p13_GENCODE/" ../0_mapping.sh > 00.mapping/s0_mapping_1.sh 
    sed "s/EXAMPLE-SN/${SN}/; s/EXAMPLE-FC/${FC}/; s/EXAMPLE-OUT/${OUTDIR}/; s/EXAMPLE-FQ/${FQDIR}/; s/GRCm39_GENCODE/GRCh38.p13_GENCODE/; s/L01/L02/" ../0_mapping.sh > 00.mapping/s0_mapping_2.sh 
    sed "s/EXAMPLE-SN/${SN}/; s/EXAMPLE-FC/${FC}/; s/EXAMPLE-OUT/${OUTDIR}/; s/EXAMPLE-FQ/${FQDIR}/; s/GRCm39_GENCODE/GRCh38.p13_GENCODE/; s/L01/L03/" ../0_mapping.sh > 00.mapping/s0_mapping_3.sh 
    sed "s/EXAMPLE-SN/${SN}/; s/EXAMPLE-FC/${FC}/; s/EXAMPLE-OUT/${OUTDIR}/; s/EXAMPLE-FQ/${FQDIR}/; s/GRCm39_GENCODE/GRCh38.p13_GENCODE/; s/L01/L04/" ../0_mapping.sh > 00.mapping/s0_mapping_4.sh 
else
    sed "s/EXAMPLE-SN/${SN}/; s/EXAMPLE-FC/${FC}/; s/EXAMPLE-OUT/${OUTDIR}/; s/EXAMPLE-FQ/${FQDIR}/" ../0_mapping.sh > 00.mapping/s0_mapping_1.sh 
    sed "s/EXAMPLE-SN/${SN}/; s/EXAMPLE-FC/${FC}/; s/EXAMPLE-OUT/${OUTDIR}/; s/EXAMPLE-FQ/${FQDIR}/; s/L01/L02/" ../0_mapping.sh > 00.mapping/s0_mapping_2.sh 
    sed "s/EXAMPLE-SN/${SN}/; s/EXAMPLE-FC/${FC}/; s/EXAMPLE-OUT/${OUTDIR}/; s/EXAMPLE-FQ/${FQDIR}/; s/L01/L03/" ../0_mapping.sh > 00.mapping/s0_mapping_3.sh 
    sed "s/EXAMPLE-SN/${SN}/; s/EXAMPLE-FC/${FC}/; s/EXAMPLE-OUT/${OUTDIR}/; s/EXAMPLE-FQ/${FQDIR}/; s/L01/L04/" ../0_mapping.sh > 00.mapping/s0_mapping_4.sh 
fi

# bcPara files
sed "s/EXAMPLE-SN/${SN}/g; s/EXAMPLE-FC/${FC}/g; s/EXAMPLE-OUT/${OUTDIR}/g; s/EXAMPLE-FQ/${FQDIR}/g" ../G400.bcPara > 00.mapping/L01.bcPara 
sed "s/EXAMPLE-SN/${SN}/g; s/EXAMPLE-FC/${FC}/g; s/EXAMPLE-OUT/${OUTDIR}/g; s/EXAMPLE-FQ/${FQDIR}/g; s/L01/L02/g" ../G400.bcPara > 00.mapping/L02.bcPara 
sed "s/EXAMPLE-SN/${SN}/g; s/EXAMPLE-FC/${FC}/g; s/EXAMPLE-OUT/${OUTDIR}/g; s/EXAMPLE-FQ/${FQDIR}/g; s/L01/L03/g" ../G400.bcPara > 00.mapping/L03.bcPara 
sed "s/EXAMPLE-SN/${SN}/g; s/EXAMPLE-FC/${FC}/g; s/EXAMPLE-OUT/${OUTDIR}/g; s/EXAMPLE-FQ/${FQDIR}/g; s/L01/L04/g" ../G400.bcPara > 00.mapping/L04.bcPara 

# merge script
sed "s/EXAMPLE-SN/${SN}/; s/EXAMPLE-FC/${FC}/; s/EXAMPLE-OUT/${OUTDIR}/" ../1_merge.sh > 01.merge/s1_merge.sh

# count script
if [[ $5 = "human" ]]; then
    sed "s/EXAMPLE-SN/${SN}/; s/EXAMPLE-FC/${FC}/; s/EXAMPLE-OUT/${OUTDIR}/; s/GRCm39_GENCODE/GRCh38.p13_GENCODE/; s/vM30/v41/" ../2_count.sh > 02.count/s2_count.sh
else
    sed "s/EXAMPLE-SN/${SN}/; s/EXAMPLE-FC/${FC}/; s/EXAMPLE-OUT/${OUTDIR}/" ../2_count.sh > 02.count/s2_count.sh
fi

# tissueCut script
sed "s/EXAMPLE-SN/${SN}/; s/EXAMPLE-OUT/${OUTDIR}/" ../4_tissueCut.sh > 04.tissuecut/s4_tissueCut.sh

# spatialCluster script
sed "s/EXAMPLE-SN/${SN}/; s/EXAMPLE-OUT/${OUTDIR}/" ../5_spatialCluster.sh > 05.spatialcluster/s5_spatialCluster.sh

# saturation script
sed "s/EXAMPLE-SN/${SN}/; s/EXAMPLE-FC/${FC}/; s/EXAMPLE-OUT/${OUTDIR}/" ../6_saturation.sh > 06.saturation/s6_saturation.sh

# report script
sed "s/EXAMPLE-SN/${SN}/; s/EXAMPLE-FC/${FC}/; s/EXAMPLE-OUT/${OUTDIR}/" ../7_report.sh > 07.report/s7_report.sh

# other scripts
sed "s/EXAMPLE-SN/${SN}/; s/EXAMPLE-OUT/${OUTDIR}/" ../CIDcount.sh > CIDcount.sh
sed "s/EXAMPLE-SN/${SN}/; s/EXAMPLE-OUT/${OUTDIR}/" ../gefTools.sh > gefTools.sh
