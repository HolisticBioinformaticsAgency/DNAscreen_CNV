#!/bin/bash
module load samtools
module load bedtools

REF=/fs04/vh83/reference/genomes/hg38/heng_li_recomended/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna
FAI=${REF}.fai
BED=/home/zlaw0001/vh83_scratch/projects/temp_dnascreen_copy/dnascreen/demultiplex_pb/workflow_HiFiTargetEnrichment/3539131_Covered_DNA_Screen_for_ONT_alignment.bed

# 1. Generate .fai if it doesn't already exist
if [ ! -f "$FAI" ]; then
    samtools faidx $REF
fi

# 2. Merge BED (in case of overlaps) and sum covered bases
COVERED=$(bedtools merge -i $BED | awk '{sum += $3 - $2} END {print sum}')
echo "Total bases covered: $COVERED"

# 3. Get total genome size from .fai next to the reference
GENOME=$(awk '{sum += $2} END {print sum}' $FAI)
echo "Total genome size: $GENOME"

# 4. Calculate percentage
echo "scale=6; $COVERED / $GENOME * 100" | bc
