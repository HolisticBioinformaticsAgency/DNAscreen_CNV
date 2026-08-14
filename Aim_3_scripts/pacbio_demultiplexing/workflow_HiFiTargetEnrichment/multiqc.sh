#!/bin/bash
#SBATCH --job-name=multiqc_pacbio
#SBATCH --output=multiqc_pacbio_%j.out
#SBATCH --error=multiqc_pacbio_%j.err
#SBATCH --time=4:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G

### SOME MINOR CHANGES NEEDED. THE LOG FILES HAVE BARCODES 10-16 THAT DON'T FIT THE 97-14 DESIGNATION
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
BASE_DIR="/home/zlaw0001/vh83_scratch/projects/temp_dnascreen_copy/dnascreen/demultiplex_pb/workflow_HiFiTargetEnrichment/batches/target_hifit"
STATS_BASE_DIR="${BASE_DIR}/pacbio_stats"
MULTIQC_OUT="${BASE_DIR}/multiqc_report"

# pbmarkdup log directory — one log file per barcode
PBMARKDUP_LOG_DIR="/fs04/scratch2/vh83/projects/temp_dnascreen_copy/dnascreen/demultiplex_pb/workflow_HiFiTargetEnrichment/batches/target_hifit/logs/pbmarkdup"

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
# PacBio HiFi Q20 threshold
# rq tag = predicted read accuracy (0-1); rq >= 0.99 = Q20 (99% accuracy)
# NOTE: if rq tag is absent in realigned_minimap2 BAMs (stripped during FASTQ
#       conversion), Q20 filtering will silently pass 0 reads.
#       The script checks for this and warns if the Q20 BAM is empty.
############################
MIN_RQ="0.99"

############################
# Barcode list
############################
BARCODES=(
  bc1001 bc1002 bc1003 bc1004 bc1005 bc1006 bc1007 bc1008
  bc1097 bc1098 bc1099 bc1100 bc1101 bc1102 bc1103 bc1104
)

# Maps sample name → actual pbmarkdup log filename stem
declare -A BC_LOG_MAP=(
  [bc1001]=bc1001  [bc1002]=bc1002  [bc1003]=bc1003  [bc1004]=bc1004
  [bc1005]=bc1005  [bc1006]=bc1006  [bc1007]=bc1007  [bc1008]=bc1008
  [bc1097]=bc1009  [bc1098]=bc1010  [bc1099]=bc1011  [bc1100]=bc1012
  [bc1101]=bc1013  [bc1102]=bc1014  [bc1103]=bc1015  [bc1104]=bc1016
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
# Helper: parse duplication rate from pbmarkdup log file
# Usage: parse_dup_rate_from_log <log_file> <stats_dir> <sample_name>
#
# Parses the TOTAL summary line in the pbmarkdup log, e.g.:
#   TOTAL    884296    812425 (91.9%)    71871 (8.1%)
#
# Extracts:
#   Total_reads_raw   = column 2 (READS)
#   Unique_reads      = column 3 (first token of UNIQUE MOLECULES)
#   Duplicate_reads   = column 4 (first token of DUPLICATE READS)
#   Duplication_rate  = decimal form of the percentage in DUPLICATE READS
############################
parse_dup_rate_from_log() {
  local LOG_FILE="$1"
  local STATS_DIR="$2"
  local SAMPLE="$3"
  local OUT_TSV="${STATS_DIR}/${SAMPLE}_dup_rate.tsv"

  echo -e "Sample\tTotal_reads_raw\tUnique_reads\tDuplicate_reads\tDuplication_rate_pct" > "${OUT_TSV}"

  if [[ ! -f "${LOG_FILE}" ]]; then
    echo "[WARNING] ${SAMPLE}: pbmarkdup log not found — dup rate will be NA: ${LOG_FILE}"
    echo -e "${SAMPLE}\tNA\tNA\tNA\tNA" >> "${OUT_TSV}"
    return
  fi

  # The TOTAL line looks like:
  # TOTAL                     884296      812425 (91.9%)       71871 (8.1%)
  local TOTAL_LINE
  TOTAL_LINE=$(grep -m1 "^TOTAL" "${LOG_FILE}" || true)

  if [[ -z "${TOTAL_LINE}" ]]; then
    echo "[WARNING] ${SAMPLE}: no TOTAL line found in pbmarkdup log — dup rate will be NA"
    echo -e "${SAMPLE}\tNA\tNA\tNA\tNA" >> "${OUT_TSV}"
    return
  fi

  # Use awk to parse: fields are whitespace-separated with parenthetical percentages
  # TOTAL  <reads>  <unique> (<pct>%)  <dups> (<dup_pct>%)
  local TOTAL_READS UNIQUE_READS DUP_READS DUP_RATE
  read -r TOTAL_READS UNIQUE_READS DUP_READS DUP_RATE <<< "$(
    echo "${TOTAL_LINE}" | awk '{
      # $1=TOTAL $2=reads $3=unique $4=(pct%) $5=dups $6=(dup_pct%)
      reads  = $2
      unique = $3
      dups   = $5
      # parse percentage from $6: "(8.1%)" → 0.0810
      pct_str = $6
      gsub(/[()%]/, "", pct_str)
      rate = (pct_str + 0) / 100
      printf "%s %s %s %.4f\n", reads, unique, dups, rate
    }'
  )"

  echo -e "${SAMPLE}\t${TOTAL_READS}\t${UNIQUE_READS}\t${DUP_READS}\t${DUP_RATE}" >> "${OUT_TSV}"
  echo "[${SAMPLE}] Duplication rate: ${DUP_RATE} (${DUP_READS} dups / ${TOTAL_READS} total)"
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
      # Compute % on-target (guard against divide-by-zero)
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
# Helper: run all stats for one realigned_minimap2 BAM
# Usage: run_stats <bam> <stats_dir> <sample_name> <pbmarkdup_log>
#
# These BAMs are treated as already deduplicated (pbmarkdup ran before alignment,
# duplicates were removed prior to FASTQ conversion → minimap2 realignment).
# Therefore NO -F 1024 filter is applied.
#
# Duplication rate    → parsed from pbmarkdup log (TOTAL line)
# Total Reads section → BAM as-is (all reads = deduped reads)
# >=Q20 Reads section → BAM filtered by rq >= 0.99 (HiFi Q20)
#                       WARNING: rq tag may be absent if stripped during
#                       FASTQ conversion. Check Q20 read count in output.
# Target Enrichment   → Picard HsMetrics on BAM, run TWICE (9genes + covered)
############################
run_stats() {
  local BAM="$1"
  local STATS_DIR="$2"
  local SAMPLE="$3"
  local PBMARKDUP_LOG="$4"

  mkdir -p "${STATS_DIR}"

  if [[ ! -f "${BAM}.bai" ]]; then
    samtools index -@ "${THREADS}" "${BAM}"
  fi

  # --- Duplication rate from pbmarkdup log ---
  parse_dup_rate_from_log "${PBMARKDUP_LOG}" "${STATS_DIR}" "${SAMPLE}"

  # --- Q20 temp BAM ---
  Q20_BAM="${STATS_DIR}/${SAMPLE}_q20.tmp.bam"

  samtools view -@ "${THREADS}" -b -e "[rq]>=${MIN_RQ}" "${BAM}" -o "${Q20_BAM}"
  samtools index -@ "${THREADS}" "${Q20_BAM}"

  Q20_COUNT=$(samtools view -c "${Q20_BAM}")
  if [[ "${Q20_COUNT}" -eq 0 ]]; then
    echo "[WARNING] ${SAMPLE}: Q20 BAM is empty — rq tag may have been stripped during FASTQ conversion."
    echo "[WARNING] ${SAMPLE}: Falling back to copying Total Reads stats for the Q20 section."
    cp "${Q20_BAM}" "${Q20_BAM}.empty_flag"
  fi

  # --- samtools stats on full BAM (Total Reads section) ---
  samtools stats -@ "${THREADS}" "${BAM}" \
    > "${STATS_DIR}/${SAMPLE}.stats.txt"

  parse_read_stats \
    "${STATS_DIR}/${SAMPLE}.stats.txt" \
    "${SAMPLE}" \
    "${STATS_DIR}/${SAMPLE}_read_summary.tsv"

  # --- samtools stats on Q20 BAM ---
  if [[ "${Q20_COUNT}" -eq 0 ]]; then
    cp "${STATS_DIR}/${SAMPLE}_read_summary.tsv" \
       "${STATS_DIR}/${SAMPLE}_read_summary_q20.tsv"
    echo "[WARNING] ${SAMPLE}: Q20 summary is a copy of Total Reads (rq tag absent)."
  else
    samtools stats -@ "${THREADS}" "${Q20_BAM}" \
      > "${STATS_DIR}/${SAMPLE}_q20.stats.txt"

    parse_read_stats \
      "${STATS_DIR}/${SAMPLE}_q20.stats.txt" \
      "${SAMPLE}" \
      "${STATS_DIR}/${SAMPLE}_read_summary_q20.tsv"
  fi

  # --- flagstat / idxstats / coverage on full BAM ---
  samtools flagstat -@ "${THREADS}" "${BAM}" \
    > "${STATS_DIR}/${SAMPLE}.flagstat.txt"

  samtools idxstats "${BAM}" \
    > "${STATS_DIR}/${SAMPLE}.idxstats.txt"

  samtools coverage "${BAM}" \
    > "${STATS_DIR}/${SAMPLE}.coverage.txt"

  # --- Picard HsMetrics: 9genes interval list ---
  run_hsmetrics "${BAM}" "${STATS_DIR}" "${SAMPLE}" "${INTERVAL_LIST_9GENES}" "9genes"

  # --- Picard HsMetrics: Covered BED interval list ---
  run_hsmetrics "${BAM}" "${STATS_DIR}" "${SAMPLE}" "${INTERVAL_LIST_COVERED}" "covered"

  # --- Cleanup temp BAMs ---
  rm -f "${Q20_BAM}" "${Q20_BAM}.bai" "${Q20_BAM}.empty_flag"
}

############################
# Per-barcode stats from realigned_minimap2 BAMs
# Duplication rate is parsed from the pre-existing pbmarkdup log files.
############################
BAM_LIST=()

for BC in "${BARCODES[@]}"; do
  echo "===================================="
  echo "Generating stats: ${BC}"
  echo "===================================="

  REALIGNED_BAM="${BASE_DIR}/${BC}/realigned_minimap2/${BC}.minimap2.GRCh38.bam"
  PBMARKDUP_LOG="${PBMARKDUP_LOG_DIR}/${BC_LOG_MAP[${BC}]}.log"
  STATS_DIR="${STATS_BASE_DIR}/${BC}"

  if [[ ! -f "${REALIGNED_BAM}" ]]; then
    echo "[WARNING] realigned BAM not found for ${BC}, skipping: ${REALIGNED_BAM}"
    continue
  fi

  run_stats "${REALIGNED_BAM}" "${STATS_DIR}" "${BC}" "${PBMARKDUP_LOG}"
  BAM_LIST+=("${REALIGNED_BAM}")

  echo "[${BC}] Stats done"
done

############################
# TOTAL (merged) stats
# Duplication rate for TOTAL is aggregated by summing raw counts across all
# per-barcode dup_rate.tsv files (avoids needing a merged markdup BAM/log).
############################
echo "===================================="
echo "Generating stats: TOTAL (all barcodes merged)"
echo "===================================="

TOTAL_STATS_DIR="${STATS_BASE_DIR}/total"
MERGED_BAM="${STATS_BASE_DIR}/total/total_merged.bam"
mkdir -p "${TOTAL_STATS_DIR}"

if [[ ${#BAM_LIST[@]} -gt 0 ]]; then
  samtools merge -f -@ "${THREADS}" "${MERGED_BAM}" "${BAM_LIST[@]}"
  samtools index -@ "${THREADS}" "${MERGED_BAM}"

  # Aggregate dup rate across all barcodes from their individual tsv files
  # (no log file exists for a merged total — we compute it from parsed counts)
  TOTAL_DUP_TSV="${TOTAL_STATS_DIR}/total_dup_rate.tsv"
  echo -e "Sample\tTotal_reads_raw\tUnique_reads\tDuplicate_reads\tDuplication_rate" > "${TOTAL_DUP_TSV}"

  awk 'BEGIN{OFS="\t"; tot=0; uniq=0; dups=0}
    FNR==1{next}  # skip header in each file
    $2 ~ /^[0-9]+$/ {
      tot  += $2
      uniq += $3
      dups += $4
    }
    END{
      rate = (tot > 0) ? dups/tot : 0
      printf "total\t%d\t%d\t%d\t%.4f\n", tot, uniq, dups, rate
    }' \
    $(find "${STATS_BASE_DIR}" -mindepth 2 -maxdepth 2 -name "*_dup_rate.tsv" ! -path "*/total/*") \
    >> "${TOTAL_DUP_TSV}"

  # Run all other stats (Q20 filter, samtools stats, Picard) on merged BAM.
  # For run_stats the 4th arg (log) is not used when dup rate TSV already exists,
  # so we pass a dummy path — parse_dup_rate_from_log will warn and write NA,
  # but we immediately overwrite with the aggregated TSV computed above.
  DUMMY_LOG="/dev/null"

  # Run remaining stats; dup rate TSV is pre-populated so we skip parse step
  # by running a slim version inline:
  STATS_DIR="${TOTAL_STATS_DIR}"
  SAMPLE="total"
  BAM="${MERGED_BAM}"

  mkdir -p "${STATS_DIR}"
  if [[ ! -f "${BAM}.bai" ]]; then
    samtools index -@ "${THREADS}" "${BAM}"
  fi

  # Q20 temp BAM
  Q20_BAM="${STATS_DIR}/${SAMPLE}_q20.tmp.bam"
  samtools view -@ "${THREADS}" -b -e "[rq]>=${MIN_RQ}" "${BAM}" -o "${Q20_BAM}"
  samtools index -@ "${THREADS}" "${Q20_BAM}"
  Q20_COUNT=$(samtools view -c "${Q20_BAM}")

  if [[ "${Q20_COUNT}" -eq 0 ]]; then
    echo "[WARNING] total: Q20 BAM is empty — rq tag may have been stripped during FASTQ conversion."
    cp "${Q20_BAM}" "${Q20_BAM}.empty_flag"
  fi

  samtools stats -@ "${THREADS}" "${BAM}" > "${STATS_DIR}/${SAMPLE}.stats.txt"
  parse_read_stats "${STATS_DIR}/${SAMPLE}.stats.txt" "${SAMPLE}" "${STATS_DIR}/${SAMPLE}_read_summary.tsv"

  if [[ "${Q20_COUNT}" -eq 0 ]]; then
    cp "${STATS_DIR}/${SAMPLE}_read_summary.tsv" "${STATS_DIR}/${SAMPLE}_read_summary_q20.tsv"
    echo "[WARNING] total: Q20 summary is a copy of Total Reads (rq tag absent)."
  else
    samtools stats -@ "${THREADS}" "${Q20_BAM}" > "${STATS_DIR}/${SAMPLE}_q20.stats.txt"
    parse_read_stats "${STATS_DIR}/${SAMPLE}_q20.stats.txt" "${SAMPLE}" "${STATS_DIR}/${SAMPLE}_read_summary_q20.tsv"
  fi

  samtools flagstat -@ "${THREADS}" "${BAM}" > "${STATS_DIR}/${SAMPLE}.flagstat.txt"
  samtools idxstats "${BAM}" > "${STATS_DIR}/${SAMPLE}.idxstats.txt"
  samtools coverage "${BAM}" > "${STATS_DIR}/${SAMPLE}.coverage.txt"

  run_hsmetrics "${BAM}" "${STATS_DIR}" "${SAMPLE}" "${INTERVAL_LIST_9GENES}" "9genes"
  run_hsmetrics "${BAM}" "${STATS_DIR}" "${SAMPLE}" "${INTERVAL_LIST_COVERED}" "covered"

  rm -f "${Q20_BAM}" "${Q20_BAM}.bai" "${Q20_BAM}.empty_flag"

  echo "[TOTAL] Stats done"
else
  echo "[WARNING] No realigned BAMs found — skipping TOTAL merge"
fi

############################
# Build MultiQC custom tables
############################
echo "===================================="
echo "Building custom MultiQC tables"
echo "===================================="

CUSTOM_DIR="${STATS_BASE_DIR}/multiqc_custom_content"
mkdir -p "${CUSTOM_DIR}"

READS_TABLE="${CUSTOM_DIR}/pb_read_summary_mqc.tsv"
Q20_TABLE="${CUSTOM_DIR}/pb_q20_read_summary_mqc.tsv"
DUP_TABLE="${CUSTOM_DIR}/pb_dup_rate_mqc.tsv"
TARGET_9GENES_TABLE="${CUSTOM_DIR}/pb_target_9genes_mqc.tsv"
TARGET_COVERED_TABLE="${CUSTOM_DIR}/pb_target_covered_mqc.tsv"
MQC_CONFIG="${STATS_BASE_DIR}/multiqc_config.yaml"

READ_SUMMARY_FILES=$(find "${STATS_BASE_DIR}" -mindepth 2 -maxdepth 2 -name "*_read_summary.tsv" ! -name "*_q20*")
Q20_SUMMARY_FILES=$(find "${STATS_BASE_DIR}" -mindepth 2 -maxdepth 2 -name "*_read_summary_q20.tsv")
DUP_RATE_FILES=$(find "${STATS_BASE_DIR}" -mindepth 2 -maxdepth 2 -name "*_dup_rate.tsv")
HSMETRICS_9GENES_FILES=$(find "${STATS_BASE_DIR}" -mindepth 2 -maxdepth 2 -name "*_hsmetrics_9genes_summary.tsv")
HSMETRICS_COVERED_FILES=$(find "${STATS_BASE_DIR}" -mindepth 2 -maxdepth 2 -name "*_hsmetrics_covered_summary.tsv")

echo -e "Sample\tTotal_reads\tTotal_bp\tN50\tMax_length_bp\tAvg_length_bp\tMedian_length_bp\tAvg_qscore" > "${READS_TABLE}"
[[ -n "${READ_SUMMARY_FILES}" ]] && awk 'FNR>1{print}' ${READ_SUMMARY_FILES} >> "${READS_TABLE}"

echo -e "Sample\tTotal_reads\tTotal_bp\tN50\tMax_length_bp\tAvg_length_bp\tMedian_length_bp\tAvg_qscore" > "${Q20_TABLE}"
[[ -n "${Q20_SUMMARY_FILES}" ]] && awk 'FNR>1{print}' ${Q20_SUMMARY_FILES} >> "${Q20_TABLE}"

echo -e "Sample\tTotal_reads_raw\tUnique_reads\tDuplicate_reads\tDuplication_rate" > "${DUP_TABLE}"
[[ -n "${DUP_RATE_FILES}" ]] && awk 'FNR>1{print}' ${DUP_RATE_FILES} >> "${DUP_TABLE}"

echo -e "Sample\tPCT_SELECTED_BASES\tMEAN_TARGET_COVERAGE\tFOLD_ENRICHMENT\tFOLD_80_BASE_PENALTY\tPCT_TARGET_BASES_10X\tPCT_TARGET_BASES_20X\tPCT_TARGET_BASES_50X\tPCT_TARGET_BASES_100X\tAT_DROPOUT\tGC_DROPOUT\tPCT_ON_TARGET" > "${TARGET_9GENES_TABLE}"
[[ -n "${HSMETRICS_9GENES_FILES}" ]] && awk 'FNR>1{print}' ${HSMETRICS_9GENES_FILES} >> "${TARGET_9GENES_TABLE}"

echo -e "Sample\tPCT_SELECTED_BASES\tMEAN_TARGET_COVERAGE\tFOLD_ENRICHMENT\tFOLD_80_BASE_PENALTY\tPCT_TARGET_BASES_10X\tPCT_TARGET_BASES_20X\tPCT_TARGET_BASES_50X\tPCT_TARGET_BASES_100X\tAT_DROPOUT\tGC_DROPOUT\tPCT_ON_TARGET" > "${TARGET_COVERED_TABLE}"
[[ -n "${HSMETRICS_COVERED_FILES}" ]] && awk 'FNR>1{print}' ${HSMETRICS_COVERED_FILES} >> "${TARGET_COVERED_TABLE}"

cat > "${MQC_CONFIG}" << 'EOF'
custom_data:
  pb_total_reads:
    file_format: "tsv"
    section_name: "Total Reads (deduped)"
    description: "Read-level summary from realigned_minimap2 BAMs, treated as deduplicated (pbmarkdup ran before FASTQ conversion)."
    plot_type: "table"
    pconfig:
      id: "pb_total_reads_table"
      namespace: "PacBio QC"
  pb_dup_rate:
    file_format: "tsv"
    section_name: "Duplication Rate"
    description: "Duplication rate parsed from pbmarkdup log files (TOTAL summary line). Duplication_rate = Duplicate_reads / Total_reads_raw."
    plot_type: "table"
    pconfig:
      id: "pb_dup_rate_table"
      namespace: "PacBio QC"
  pb_q20_reads:
    file_format: "tsv"
    section_name: ">= Q20 HiFi Reads (deduped)"
    description: "Read-level summary filtered to rq >= 0.99 (HiFi Q20). Note: if rq tag was stripped during FASTQ conversion, values will match Total Reads."
    plot_type: "table"
    pconfig:
      id: "pb_q20_reads_table"
      namespace: "PacBio QC"
  pb_target_9genes:
    file_format: "tsv"
    section_name: "Target Enrichment — 9genes (deduped)"
    description: "Picard CollectHsMetrics. PCT_ON_TARGET = ON_TARGET_BASES / PF_BASES_ALIGNED × 100."
    plot_type: "table"
    pconfig:
      id: "pb_target_9genes_table"
      namespace: "PacBio QC"
      col_config:
        PCT_ON_TARGET:
          title: "% On-Target"
          description: "Percentage of aligned bases that fall on a target region (ON_TARGET_BASES / PF_BASES_ALIGNED × 100)"
          min: 0
          max: 100
          suffix: "%"
  pb_target_covered:
    file_format: "tsv"
    section_name: "Target Enrichment — Covered BED (deduped)"
    description: "Percentage of aligned bases that fall on a target region (ON_TARGET_BASES / PF_BASES_ALIGNED × 100)"
    plot_type: "table"
    pconfig:
      id: "pb_target_covered_table"
      namespace: "PacBio QC"
      col_config:
        PCT_ON_TARGET:
          title: "% On-Target"
          description: "Percentage of aligned bases that fall on a target region (ON_TARGET_BASES / PF_BASES_ALIGNED × 100)"
          min: 0
          max: 100
          suffix: "%"
sp:
  pb_total_reads:
    fn: "pb_read_summary_mqc.tsv"
  pb_dup_rate:
    fn: "pb_dup_rate_mqc.tsv"
  pb_q20_reads:
    fn: "pb_q20_read_summary_mqc.tsv"
  pb_target_9genes:
    fn: "pb_target_9genes_mqc.tsv"
  pb_target_covered:
    fn: "pb_target_covered_mqc.tsv"
custom_content:
  order:
    - pb_total_reads
    - pb_dup_rate
    - pb_q20_reads
    - pb_target_9genes
    - pb_target_covered
EOF

############################
# Run MultiQC
############################
echo "===================================="
echo "Running MultiQC"
echo "===================================="

mkdir -p "${MULTIQC_OUT}"

multiqc \
  --config "${MQC_CONFIG}" \
  --outdir "${MULTIQC_OUT}" \
  --filename "pacbio_multiqc_report" \
  --title "PacBio HiFi Target Enrichment QC" \
  --ignore "*/multiqc_report*" \
  --ignore "*/multiqc_data*" \
  "${STATS_BASE_DIR}"

echo "MultiQC report written to: ${MULTIQC_OUT}/pacbio_multiqc_report.html 🎉"
