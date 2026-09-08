#!/bin/bash
#SBATCH --job-name=multiqc_ont
#SBATCH --output=multiqc_ont_%j.out
#SBATCH --error=multiqc_ont_%j.err
#SBATCH --time=4:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G

### CHANGE DUP RATE TO 2 DECIMAL POINTS

set -euo pipefail

############################
# Load modules
############################
module load samtools
module load picard
module load multiqc

############################
# SLURM resources
############################
THREADS=${SLURM_CPUS_PER_TASK}

############################
# Paths
############################
BASE_DIR="/home/zlaw0001/vh83_scratch/projects/temp_dnascreen_copy/dnascreen/dorado_demultiplex_142pods_barcode_both_ends/bam_stats_dedup"
MULTIQC_OUT="${BASE_DIR}/multiqc_report_bams_dedup_cleaned"

# Source BAM directories
MARKDUP_DIR="/home/zlaw0001/vh83_scratch/projects/temp_dnascreen_copy/dnascreen/dorado_demultiplex_142pods_barcode_both_ends/bams_markdup"

############################
# Reference & target files
#
# Two interval lists — both run for every BAM:
#   INTERVAL_LIST_9GENES  → "9genes"      (DNA Screen ONT panel)
#   INTERVAL_LIST_COVERED → "Covered BED" (full covered BED)
############################
REF_FASTA="/fs04/vh83/reference/genomes/hg38/heng_li_recomended/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna"
INTERVAL_LIST_9GENES="/fs04/scratch2/vh83/projects/temp_dnascreen_copy/dnascreen/dorado_demultiplex_142pods_barcode_both_ends/3539131_Covered_DNA_Screen_for_ONT_alignment.interval_list"
INTERVAL_LIST_COVERED="/fs04/scratch2/vh83/projects/temp_dnascreen_copy/dnascreen/dorado_demultiplex_142pods_barcode_both_ends/3539131_Covered.interval_list"

############################
# Barcode list
############################
BARCODES=(
  barcode1001 barcode1002 barcode1003 barcode1004 barcode1005 barcode1006 barcode1007 barcode1008
  barcode1097 barcode1098 barcode1099 barcode1100 barcode1101 barcode1102 barcode1103 barcode1104
)

############################
# Validate required files
############################
for F in "${REF_FASTA}" "${INTERVAL_LIST_9GENES}" "${INTERVAL_LIST_COVERED}"; do
  if [[ ! -f "${F}" ]]; then
    echo "[ERROR] Required file not found: ${F}"
    exit 1
  fi
done

############################
# Helper: compute duplication rate from a markdup BAM (before -F 1024 filter)
# Duplication_rate = duplicates / total_reads (from samtools flagstat on raw BAM)
# Writes: <sample>_dup_rate.tsv
# Usage: parse_dup_rate <raw_bam> <stats_dir> <sample_name>
############################
parse_dup_rate() {
  local RAW_BAM="$1"
  local STATS_DIR="$2"
  local SAMPLE="$3"
  local OUT_TSV="${STATS_DIR}/${SAMPLE}_dup_rate.tsv"

  local FLAGSTAT_TMP="${STATS_DIR}/${SAMPLE}_raw.flagstat.txt"

  # Run flagstat on the raw markdup BAM (duplicates still flagged)
  samtools flagstat -@ "${THREADS}" "${RAW_BAM}" > "${FLAGSTAT_TMP}"

  # Extract total reads and duplicate-flagged reads
  local TOTAL DUPS DUP_RATE
  TOTAL=$(awk '/in total/{print $1; exit}' "${FLAGSTAT_TMP}")
  DUPS=$(awk '/duplicates/{print $1; exit}' "${FLAGSTAT_TMP}")

  if [[ "${TOTAL}" -gt 0 ]]; then
    DUP_RATE=$(awk -v d="${DUPS}" -v t="${TOTAL}" 'BEGIN{printf "%.4f", d/t}')
  else
    DUP_RATE="NA"
  fi

  echo -e "Sample\tTotal_reads_raw\tDuplicate_reads\tDuplication_rate" > "${OUT_TSV}"
  echo -e "${SAMPLE}\t${TOTAL}\t${DUPS}\t${DUP_RATE}" >> "${OUT_TSV}"

  # Clean up temp flagstat
  rm -f "${FLAGSTAT_TMP}"
}

############################
# Helper: parse read stats from a samtools stats file
# Usage: parse_read_stats <stats_file> <sample_name> <output_tsv>
############################
parse_read_stats() {
  local STATS_FILE="$1"
  local SAMPLE="$2"
  local OUT_TSV="$3"

  local TOTAL_READS TOTAL_BP MAX_LEN AVG_LEN AVG_QUAL N50 MEDIAN_LEN

  TOTAL_READS=$(grep "^SN" "${STATS_FILE}" | awk -F'\t' '/raw total sequences/{print $3}')
  TOTAL_BP=$(grep    "^SN" "${STATS_FILE}" | awk -F'\t' '/total length/{print $3}'        | head -1)
  MAX_LEN=$(grep     "^SN" "${STATS_FILE}" | awk -F'\t' '/maximum length/{print $3}')
  AVG_LEN=$(grep     "^SN" "${STATS_FILE}" | awk -F'\t' '/average length/{print $3}')
  AVG_QUAL=$(grep    "^SN" "${STATS_FILE}" | awk -F'\t' '/average quality/{print $3}')

  N50=$(awk '
    /^RL/ { rl[$2]+=$3; total+=$2*$3 }
    END {
      half=total/2; run=0;
      n=asorti(rl,keys)
      for(i=n;i>=1;i--) {
        run+=keys[i]*rl[keys[i]]
        if(run>=half){ print keys[i]; break }
      }
    }' "${STATS_FILE}")

  MEDIAN_LEN=$(awk '
    /^RL/ { rl[$2]+=$3; total+=$3 }
    END {
      target=(total+1)/2; run=0;
      n=asorti(rl,keys)
      for(i=1;i<=n;i++) {
        run+=rl[keys[i]]
        if(run>=target){ print keys[i]; break }
      }
    }' "${STATS_FILE}")

  echo -e "Sample\tTotal_reads\tTotal_bp\tN50\tMax_length_bp\tAvg_length_bp\tMedian_length_bp\tAvg_qscore" > "${OUT_TSV}"
  echo -e "${SAMPLE}\t${TOTAL_READS}\t${TOTAL_BP}\t${N50}\t${MAX_LEN}\t${AVG_LEN}\t${MEDIAN_LEN}\t${AVG_QUAL}" >> "${OUT_TSV}"
}

############################
# Helper: parse Picard HsMetrics into simple TSV
# Usage: parse_hsmetrics <hsmetrics_txt> <sample_name> <output_tsv>
############################
parse_hsmetrics() {
  local HSM_FILE="$1"
  local SAMPLE="$2"
  local OUT_TSV="$3"

  awk -v sample="${SAMPLE}" 'BEGIN{FS="\t"; OFS="\t"}
    /^BAIT_SET/ {
      for (i=1; i<=NF; i++) idx[$i]=i;
      next
    }
    /^#/ { next }
    NF > 1 && $1 != "BAIT_SET" {
      on_target  = $idx["ON_TARGET_BASES"]
      pf_aligned = $idx["PF_BASES_ALIGNED"]
      pct_on_target = (pf_aligned > 0) ? (on_target / pf_aligned * 100) : 0

      print "Sample","PCT_SELECTED_BASES","MEAN_TARGET_COVERAGE","FOLD_ENRICHMENT","FOLD_80_BASE_PENALTY","PCT_TARGET_BASES_10X","PCT_TARGET_BASES_20X","PCT_TARGET_BASES_50X","PCT_TARGET_BASES_100X","AT_DROPOUT","GC_DROPOUT","PCT_ON_TARGET";
      printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%.2f\n",
        sample,
        $idx["PCT_SELECTED_BASES"],
        $idx["MEAN_TARGET_COVERAGE"],
        $idx["FOLD_ENRICHMENT"],
        $idx["FOLD_80_BASE_PENALTY"],
        $idx["PCT_TARGET_BASES_10X"],
        $idx["PCT_TARGET_BASES_20X"],
        $idx["PCT_TARGET_BASES_50X"],
        $idx["PCT_TARGET_BASES_100X"],
        $idx["AT_DROPOUT"],
        $idx["GC_DROPOUT"],
        pct_on_target;
      exit
    }' "${HSM_FILE}" > "${OUT_TSV}"
}

############################
# Helper: run Picard HsMetrics for one interval list
# Usage: run_hsmetrics <bam> <stats_dir> <sample_name> <interval_list> <label>
############################
run_hsmetrics() {
  local BAM="$1"
  local STATS_DIR="$2"
  local SAMPLE="$3"
  local INTERVAL_LIST="$4"
  local LABEL="$5"

  picard CollectHsMetrics \
    I="${BAM}" \
    O="${STATS_DIR}/${SAMPLE}.HsMetrics_${LABEL}.txt" \
    R="${REF_FASTA}" \
    BAIT_INTERVALS="${INTERVAL_LIST}" \
    TARGET_INTERVALS="${INTERVAL_LIST}" \
    COVERAGE_CAP=5000 \
    2> "${STATS_DIR}/${SAMPLE}.HsMetrics_${LABEL}.log"

  parse_hsmetrics \
    "${STATS_DIR}/${SAMPLE}.HsMetrics_${LABEL}.txt" \
    "${SAMPLE}" \
    "${STATS_DIR}/${SAMPLE}_hsmetrics_${LABEL}_summary.tsv"
}

############################
# Helper: run all stats for one markdup BAM
# Usage: run_stats <markdup_bam> <stats_dir> <sample_name>
#
# Duplication rate    → computed from raw markdup BAM (before -F 1024)
# Total Reads section → markdup BAM with -F 1024 (exclude duplicates)
# >=Q10 Reads section → same BAM with -F 1024 AND qs:f >= 10
# Target Enrichment   → Picard on dedup-only BAM, run TWICE (9genes + covered)
############################
run_stats() {
  local BAM="$1"
  local STATS_DIR="$2"
  local SAMPLE="$3"

  mkdir -p "${STATS_DIR}"

  if [[ ! -f "${BAM}.bai" ]]; then
    samtools index -@ "${THREADS}" "${BAM}"
  fi

  # --- Duplication rate from raw markdup BAM (before filtering) ---
  parse_dup_rate "${BAM}" "${STATS_DIR}" "${SAMPLE}"

  # --- Temp BAMs ---
  DEDUP_BAM="${STATS_DIR}/${SAMPLE}_dedup_only.tmp.bam"
  Q10_BAM="${STATS_DIR}/${SAMPLE}_dedup_q10.tmp.bam"

  # Dedup-only: exclude duplicate flag 1024, no Q filter
  samtools view -@ "${THREADS}" -b -F 1024 "${BAM}" -o "${DEDUP_BAM}"
  samtools index -@ "${THREADS}" "${DEDUP_BAM}"

  # Q10+: exclude duplicates AND filter qs:f >= 10
  samtools view -h -@ "${THREADS}" -F 1024 "${BAM}" \
    | awk 'BEGIN{OFS="\t"}
        /^@/ { print; next }
        {
          for (i=12; i<=NF; i++) {
            if ($i ~ /^qs:f:/) {
              split($i,a,":")
              if (a[3] >= 10) print
              break
            }
          }
        }' \
    | samtools view -@ "${THREADS}" -b -o "${Q10_BAM}"
  samtools index -@ "${THREADS}" "${Q10_BAM}"

  # --- samtools stats: Total Reads (deduped) ---
  samtools stats -@ "${THREADS}" "${DEDUP_BAM}" \
    > "${STATS_DIR}/${SAMPLE}.stats.txt"

  parse_read_stats \
    "${STATS_DIR}/${SAMPLE}.stats.txt" \
    "${SAMPLE}" \
    "${STATS_DIR}/${SAMPLE}_read_summary.tsv"

  # --- samtools stats: >=Q10 Reads (deduped) ---
  samtools stats -@ "${THREADS}" "${Q10_BAM}" \
    > "${STATS_DIR}/${SAMPLE}_q10.stats.txt"

  parse_read_stats \
    "${STATS_DIR}/${SAMPLE}_q10.stats.txt" \
    "${SAMPLE}" \
    "${STATS_DIR}/${SAMPLE}_read_summary_q10.tsv"

  # --- flagstat / idxstats / coverage ---
  samtools flagstat -@ "${THREADS}" "${DEDUP_BAM}" \
    > "${STATS_DIR}/${SAMPLE}.flagstat.txt"

  samtools idxstats "${DEDUP_BAM}" \
    > "${STATS_DIR}/${SAMPLE}.idxstats.txt"

  samtools coverage "${DEDUP_BAM}" \
    > "${STATS_DIR}/${SAMPLE}.coverage.txt"

  # --- Picard HsMetrics: 9genes interval list ---
  run_hsmetrics "${DEDUP_BAM}" "${STATS_DIR}" "${SAMPLE}" "${INTERVAL_LIST_9GENES}" "9genes"

  # --- Picard HsMetrics: Covered BED interval list ---
  run_hsmetrics "${DEDUP_BAM}" "${STATS_DIR}" "${SAMPLE}" "${INTERVAL_LIST_COVERED}" "covered"

  # --- Cleanup temp BAMs ---
  rm -f "${DEDUP_BAM}" "${DEDUP_BAM}.bai" "${Q10_BAM}" "${Q10_BAM}.bai"
}

############################
# Per-barcode stats from markdup BAMs
############################
BAM_LIST=()

for barcode in "${BARCODES[@]}"; do
  echo "===================================="
  echo "Generating stats: ${barcode}"
  echo "===================================="

  MARKDUP_BAM="${MARKDUP_DIR}/markdup_${barcode}.bam"
  STATS_DIR="${BASE_DIR}/${barcode}"

  if [[ ! -f "${MARKDUP_BAM}" ]]; then
    echo "[WARNING] markdup BAM not found for ${barcode}, skipping: ${MARKDUP_BAM}"
    continue
  fi

  run_stats "${MARKDUP_BAM}" "${STATS_DIR}" "${barcode}"
  BAM_LIST+=("${MARKDUP_BAM}")

  echo "[${barcode}] Stats done"
done

############################
# TOTAL (merged) stats
############################
echo "===================================="
echo "Generating stats: TOTAL (all barcodes merged)"
echo "===================================="

TOTAL_STATS_DIR="${BASE_DIR}/total"
MERGED_BAM="${BASE_DIR}/total/total_merged_markdup.bam"
mkdir -p "${TOTAL_STATS_DIR}"

if [[ ${#BAM_LIST[@]} -gt 0 ]]; then
  samtools merge -f -@ "${THREADS}" "${MERGED_BAM}" "${BAM_LIST[@]}"
  samtools index -@ "${THREADS}" "${MERGED_BAM}"
  run_stats "${MERGED_BAM}" "${TOTAL_STATS_DIR}" "total"
  echo "[TOTAL] Stats done"
else
  echo "[WARNING] No markdup BAMs found — skipping TOTAL merge"
fi

############################
# Build MultiQC custom tables
############################
echo "===================================="
echo "Building custom MultiQC tables"
echo "===================================="

CUSTOM_DIR="${BASE_DIR}/multiqc_custom_content"
mkdir -p "${CUSTOM_DIR}"

READS_TABLE="${CUSTOM_DIR}/ont_read_summary_mqc.tsv"
Q10_TABLE="${CUSTOM_DIR}/ont_q10_read_summary_mqc.tsv"
DUP_TABLE="${CUSTOM_DIR}/ont_dup_rate_mqc.tsv"
TARGET_9GENES_TABLE="${CUSTOM_DIR}/ont_target_9genes_mqc.tsv"
TARGET_COVERED_TABLE="${CUSTOM_DIR}/ont_target_covered_mqc.tsv"
MQC_CONFIG="${BASE_DIR}/multiqc_config.yaml"

READ_SUMMARY_FILES=$(find "${BASE_DIR}" -mindepth 2 -maxdepth 2 -name "*_read_summary.tsv" ! -name "*_q10*")
Q10_SUMMARY_FILES=$(find "${BASE_DIR}" -mindepth 2 -maxdepth 2 -name "*_read_summary_q10.tsv")
DUP_RATE_FILES=$(find "${BASE_DIR}" -mindepth 2 -maxdepth 2 -name "*_dup_rate.tsv")
HSMETRICS_9GENES_FILES=$(find "${BASE_DIR}" -mindepth 2 -maxdepth 2 -name "*_hsmetrics_9genes_summary.tsv")
HSMETRICS_COVERED_FILES=$(find "${BASE_DIR}" -mindepth 2 -maxdepth 2 -name "*_hsmetrics_covered_summary.tsv")

echo -e "Sample\tTotal_reads\tTotal_bp\tN50\tMax_length_bp\tAvg_length_bp\tMedian_length_bp\tAvg_qscore" > "${READS_TABLE}"
[[ -n "${READ_SUMMARY_FILES}" ]] && awk 'FNR>1{print}' ${READ_SUMMARY_FILES} >> "${READS_TABLE}"

echo -e "Sample\tTotal_reads\tTotal_bp\tN50\tMax_length_bp\tAvg_length_bp\tMedian_length_bp\tAvg_qscore" > "${Q10_TABLE}"
[[ -n "${Q10_SUMMARY_FILES}" ]] && awk 'FNR>1{print}' ${Q10_SUMMARY_FILES} >> "${Q10_TABLE}"

echo -e "Sample\tTotal_reads_raw\tDuplicate_reads\tDuplication_rate" > "${DUP_TABLE}"
[[ -n "${DUP_RATE_FILES}" ]] && awk 'FNR>1{print}' ${DUP_RATE_FILES} >> "${DUP_TABLE}"

echo -e "Sample\tPCT_SELECTED_BASES\tMEAN_TARGET_COVERAGE\tFOLD_ENRICHMENT\tFOLD_80_BASE_PENALTY\tPCT_TARGET_BASES_10X\tPCT_TARGET_BASES_20X\tPCT_TARGET_BASES_50X\tPCT_TARGET_BASES_100X\tAT_DROPOUT\tGC_DROPOUT\tPCT_ON_TARGET" > "${TARGET_9GENES_TABLE}"
[[ -n "${HSMETRICS_9GENES_FILES}" ]] && awk 'FNR>1{print}' ${HSMETRICS_9GENES_FILES} >> "${TARGET_9GENES_TABLE}"

echo -e "Sample\tPCT_SELECTED_BASES\tMEAN_TARGET_COVERAGE\tFOLD_ENRICHMENT\tFOLD_80_BASE_PENALTY\tPCT_TARGET_BASES_10X\tPCT_TARGET_BASES_20X\tPCT_TARGET_BASES_50X\tPCT_TARGET_BASES_100X\tAT_DROPOUT\tGC_DROPOUT\tPCT_ON_TARGET" > "${TARGET_COVERED_TABLE}"
[[ -n "${HSMETRICS_COVERED_FILES}" ]] && awk 'FNR>1{print}' ${HSMETRICS_COVERED_FILES} >> "${TARGET_COVERED_TABLE}"

cat > "${MQC_CONFIG}" << 'EOF'
custom_data:
  ont_total_reads:
    file_format: "tsv"
    section_name: "Total Reads (deduped)"
    description: "Read-level summary from deduped reads (duplicates excluded via -F 1024, no Q filter)."
    plot_type: "table"
    pconfig:
      id: "ont_total_reads_table"
      namespace: "ONT QC"
  ont_dup_rate:
    file_format: "tsv"
    section_name: "Duplication Rate"
    description: "Duplication rate computed from the raw markdup BAM before duplicate removal. Duplication_rate = Duplicate_reads / Total_reads_raw (samtools flagstat 0x400 flag)."
    plot_type: "table"
    pconfig:
      id: "ont_dup_rate_table"
      namespace: "ONT QC"
  ont_q10_reads:
    file_format: "tsv"
    section_name: ">= 10 Mean_Q-Score Reads (deduped)"
    description: "Read-level summary from deduped reads additionally filtered to qs:f >= 10."
    plot_type: "table"
    pconfig:
      id: "ont_q10_reads_table"
      namespace: "ONT QC"
  ont_target_9genes:
    file_format: "tsv"
    section_name: "Target Enrichment — 9genes (deduped)"
    description: "Picard CollectHsMetrics using the DNA Screen ONT panel interval list. PCT_ON_TARGET = ON_TARGET_BASES / PF_BASES_ALIGNED × 100."
    plot_type: "table"
    pconfig:
      id: "ont_target_9genes_table"
      namespace: "ONT QC"
      col_config:
        PCT_ON_TARGET:
          title: "% On-Target"
          description: "Percentage of aligned bases on a target region (ON_TARGET_BASES / PF_BASES_ALIGNED × 100)"
          min: 0
          max: 100
          suffix: "%"
  ont_target_covered:
    file_format: "tsv"
    section_name: "Target Enrichment — Covered BED (deduped)"
    description: "Picard CollectHsMetrics using the full covered BED interval list. PCT_ON_TARGET = ON_TARGET_BASES / PF_BASES_ALIGNED × 100."
    plot_type: "table"
    pconfig:
      id: "ont_target_covered_table"
      namespace: "ONT QC"
      col_config:
        PCT_ON_TARGET:
          title: "% On-Target"
          description: "Percentage of aligned bases on a target region (ON_TARGET_BASES / PF_BASES_ALIGNED × 100)"
          min: 0
          max: 100
          suffix: "%"
sp:
  ont_total_reads:
    fn: "ont_read_summary_mqc.tsv"
  ont_dup_rate:
    fn: "ont_dup_rate_mqc.tsv"
  ont_q10_reads:
    fn: "ont_q10_read_summary_mqc.tsv"
  ont_target_9genes:
    fn: "ont_target_9genes_mqc.tsv"
  ont_target_covered:
    fn: "ont_target_covered_mqc.tsv"
custom_content:
  order:
    - ont_total_reads
    - ont_dup_rate
    - ont_q10_reads
    - ont_target_9genes
    - ont_target_covered
EOF

############################
# Run MultiQC across all
############################
echo "===================================="
echo "Running MultiQC"
echo "===================================="

mkdir -p "${MULTIQC_OUT}"

multiqc \
  --config "${MQC_CONFIG}" \
  --outdir "${MULTIQC_OUT}" \
  --filename "ONT_multiqc_report" \
  --title "ONT RUN QC" \
  --ignore "*/multiqc_report*" \
  --ignore "*/multiqc_data*" \
  "${BASE_DIR}"

echo "MultiQC report written to: ${MULTIQC_OUT}/ONT_multiqc_report.html 🎉"
