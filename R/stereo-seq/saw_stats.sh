#!/bin/bash

TOTAL_READS="$(cat *barcodeMap.stat | grep total_reads | awk '{print $2}' | paste -sd+ - | bc)"

MAPPED_CID_READS="$(cat *barcodeMap.stat | grep mapped_reads | awk '{print $2}' | paste -sd+ - | bc)"

DISCARDED_MID_READS="$(cat *barcodeMap.stat | grep umi_filter_reads | awk '{print $2}' | paste -sd+ - | bc)"

VALID_CID_READS="$(expr $MAPPED_CID_READS - $DISCARDED_MID_READS)"

INVALID_READS="$(expr $TOTAL_READS - $MAPPED_CID_READS)"

VALID_RATIO="$(expr $VALID_CID_READS / $TOTAL_READS)"

CLEAN_READS="$(cat *.final.out | grep 'Number of input reads' | awk -F '\t' '{print $2}' | paste -sd+ - | bc)"

TOO_SHORT_READS="$(expr $VALID_CID_READS - $CLEAN_READS)"

UNIQUE_MAPPING_READS="$(sed '3q;d' *.bam.summary.stat | awk '{print $2}')"

UNMAPPING_READS="$(cat *.final.out | grep 'Number of reads unmapped:' | awk -F '\t' '{print $2}' | paste -sd+ - | bc)"

MULTIPLE_MAPPING_READS="$(expr $CLEAN_READS - $UNIQUE_MAPPING_READS - $UNMAPPING_READS)"
#MULTIPLE_MAPPING_READS="$(cat *.final.out | grep 'Number of reads mapped to multiple loci' | awk -F '\t' '{print $2}' | paste -sd+ - | bc)"

READS_MAPPED_TO_GENOME="$(expr $UNIQUE_MAPPING_READS + $MULTIPLE_MAPPING_READS)"
ANNOTATION_SUCCEED="$(sed '3q;d' *.bam.summary.stat | awk '{print $3}')"
ANNOTATION_FAILED="$(expr $UNIQUE_MAPPING_READS - $ANNOTATION_SUCCEED)"

UNIQUE_READS="$(sed '3q;d' *.bam.summary.stat | awk '{print $4}')"
DUPLICATED_READS="$(expr $ANNOTATION_SUCCEED - $UNIQUE_READS)"

DUPLICATION_RATE="$(sed '3q;d' *.bam.summary.stat | awk '{print $7}')"


EXONIC_READS="$(sed '6q;d' *.bam.summary.stat | awk '{print $3}')"
INTRONIC_READS="$(sed '6q;d' *.bam.summary.stat | awk '{print $4}')"
INTERGENIC_READS="$(sed '6q;d' *.bam.summary.stat | awk '{print $5}')"

echo "--- Key Metrices ---"
#echo "---BARCODE MAPPING---"
echo "Total Reads: $TOTAL_READS"
echo "- Invalid CID Reads: $INVALID_READS"
echo "- Valid CID Reads: $VALID_CID_READS" # ($((VALID_CID_READS / TOTAL_READS)))"
echo " - Non-Relevant Short Reads: $TOO_SHORT_READS"
echo " - Clean Reads: $CLEAN_READS"
#echo " "
#echo "---GENOME MAPPING---"
echo "  - Unmapped Reads: $UNMAPPING_READS"
echo "  - Reads Mapped To Genome: $READS_MAPPED_TO_GENOME"
echo "   - Multi-Mapped Reads: $MULTIPLE_MAPPING_READS"
echo "   - Uniquely Mapped Reads: $UNIQUE_MAPPING_READS"
#echo " "
#echo "---ANNOTATION---"
echo "    - Unannotated Reads: $ANNOTATION_FAILED"
echo "    - Annotated Reads: $ANNOTATION_SUCCEED"
#echo "     - Exonic: $EXONIC_READS"
#echo "     - Intronic: $INTRONIC_READS"
#echo "     - Intergenic: $INTERGENIC_READS"
# annotated = exonic + intronic; annotated + intergenic = uniquely mapped
echo "     - Unique Reads: $UNIQUE_READS"
echo "     - Duplicated Reads: $DUPLICATED_READS"
echo " "
echo "Duplication Rate: ${DUPLICATION_RATE}%"
echo " "
echo "------"
echo "Total Reads: Total number of sequenced reads"
echo "Valid CID Reads: Number of reads with CIDs matching the mask file and with MIDs passing QC"
#echo "Invalid CID Reads: Number of reads with CIDs that cannot be matched with the mask file"
echo "Clean Reads: Number of Valid CID Reads that have passed QC"
echo "Uniquely Mapped Reads: Number of reads that mapped uniquely to the reference genome"
echo "Annotation Reads: Number of reads that are aligned to transcripts of at least one gene"
echo "Unique Reads: Number of reads which the same DNA fragment occurs and is sequenced only once"


echo "------"
echo "tissue bin stats"
grep 'binSize\|gene_type\|Umi\|^$' tissuecut.stat


Rscript ~/spatial/analysis-test/saw_sunburst.R $TOTAL_READS $INVALID_READS $VALID_CID_READS $TOO_SHORT_READS $CLEAN_READS $UNMAPPING_READS $READS_MAPPED_TO_GENOME $MULTIPLE_MAPPING_READS $UNIQUE_MAPPING_READS $ANNOTATION_FAILED $ANNOTATION_SUCCEED $UNIQUE_READS $DUPLICATED_READS
