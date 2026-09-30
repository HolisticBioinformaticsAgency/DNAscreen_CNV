#!/usr/bin/env nextflow

// ======================================================================
// Optimised CNV detection framework (per sequencing run)
//
//   1. Identify and remove outlier BAMs  (correlation -> MDS -> DBSCAN)
//   2. Call CNVs with DECoN and clearCNV on the retained BAMs
//   3. Preliminary filtering             (DECoN BF; read depth for single-exon calls)
//   4. Prioritisation matrix             (DECoN + clearCNV reciprocal overlap)
//   5. Read ratio + VAF assessment       (single-caller calls)
//   6. Plots + cross-run tracking tables
//
// Input samplesheet (CSV with header): run_id,bam_dir,vcf_dir
// ======================================================================

nextflow.enable.dsl = 2

def utils = "${projectDir}/bin/framework_utils.R"

// ----------------------------------------------------------------------
// Stage 1: read counts + outlier BAMs
// ----------------------------------------------------------------------

process COUNT_READS {

    tag "$runID"
    label 'r_small'

    publishDir "${params.outdir}/runs/${runID}/01_counts", mode: 'copy'

    input:
    tuple val(runID), val(bam_dir), val(vcf_dir)

    output:
    tuple val(runID), path("${runID}_counts.RData")

    script:
    """
    module load R

    DECON_FOLDER=\$(yq -r '.deconFolder' ${params.decon_params_nf})
    FASTA=\$(yq -r '.fasta_file' ${params.decon_params_nf})

    ls "${bam_dir}"/*.bam > all_bams.txt
    if [ ! -s all_bams.txt ]; then
        echo "No BAM files found in ${bam_dir}" >&2
        exit 1
    fi

    Rscript "\$DECON_FOLDER/ReadInBams.R" \\
        --bams all_bams.txt \\
        --bed ${params.framework_bed} \\
        --fasta "\$FASTA" \\
        --out ${runID}_counts
    """
}

process DETECT_OUTLIERS {

    tag "$runID"
    label 'r_small'

    publishDir "${params.outdir}/runs/${runID}/01_outliers", mode: 'copy', pattern: "*.{tsv,txt,rds}"
    publishDir "${params.outdir}/plots/outliers",            mode: 'copy', pattern: "*.{png,html}"

    input:
    tuple val(runID), path(counts)

    output:
    tuple val(runID), path("${runID}_outliers.tsv"),      emit: outliers
    tuple val(runID), path("${runID}_retained_bams.txt"), emit: retained
    path "${runID}_mds_coords.tsv",                       emit: mds
    path "*.{png,html,rds}",                              optional: true
    path "*_files",                                       optional: true

    script:
    """
    module load R pandoc
    detect_outliers.R --counts ${counts} --run ${runID} --minpts ${params.dbscan_minpts} --utils ${utils}
    """
}

// ----------------------------------------------------------------------
// Stage 2: CNV calling on retained BAMs
// A caller failure does not stop the run: an empty call set is emitted and
// the failure is recorded in <run>_<caller>_caller_status.tsv (-> run_status.tsv).
// ----------------------------------------------------------------------

process RUN_DECON {

    tag "$runID"

    publishDir "${params.outdir}/runs/${runID}/02_decon", mode: 'copy'

    input:
    tuple val(runID), path(retained_bams), path(counts)

    output:
    tuple val(runID), path("${runID}_decon_cnv_calls_final.tsv"), emit: calls
    path "${runID}_decon_caller_status.tsv",                     emit: status
    path "${runID}_decon.log"

    script:
    def reuse_counts = file(params.caller_bed).toAbsolutePath() == file(params.framework_bed).toAbsolutePath()
    """
    module load R

    CONFIG=${params.decon_params_nf}
    DECON_FOLDER=\$(yq -r '.deconFolder' \$CONFIG)
    FASTA=\$(yq -r '.fasta_file' \$CONFIG)
    MINCORR=\$(yq -r '.mincorr' \$CONFIG)
    MINCOV=\$(yq -r '.mincov' \$CONFIG)
    TRANSPROB=\$(yq -r '.transProb' \$CONFIG)

    OUTDIR="${runID}"
    mkdir -p "\$OUTDIR"
    STATUS="ok"
    N_BAMS=\$(grep -c . ${retained_bams} || true)

    # (set -e is suspended inside `if !`, so every step checks its own exit code)
    run_decon() {
        if [ "${reuse_counts}" = "true" ]; then
            echo "Reusing stage-1 counts (caller BED == framework BED)"
            Rscript ${projectDir}/bin/subset_counts.R ${counts} ${retained_bams} "\$OUTDIR/output.bams.RData" || return 1
        else
            Rscript "\$DECON_FOLDER/ReadInBams.R" --bams ${retained_bams} --bed ${params.caller_bed} \\
                --fasta "\$FASTA" --out "\$OUTDIR/output.bams" || return 1
        fi
        Rscript "\$DECON_FOLDER/IdentifyFailures.R" --RData "\$OUTDIR/output.bams.RData" \\
            --mincorr "\$MINCORR" --mincov "\$MINCOV" --out "\$OUTDIR/failures" || return 1
        Rscript "\$DECON_FOLDER/makeCNVcalls.R" --RData "\$OUTDIR/output.bams.RData" \\
            --transProb "\$TRANSPROB" --plot None --out "\$OUTDIR" \\
            --failures "\$OUTDIR/failures_Failures.txt" || return 1
    }

    if ! run_decon > ${runID}_decon.log 2>&1; then
        STATUS="failed"
    fi

    CNV_FILE=\$(find . -name "*_all.txt" -type f | head -n 1)
    if [ -n "\$CNV_FILE" ]; then
        cp "\$CNV_FILE" ${runID}_decon_cnv_calls_final.tsv
    else
        [ "\$STATUS" = "ok" ] && STATUS="no_calls_file"
        : > ${runID}_decon_cnv_calls_final.tsv
    fi

    printf "run\\tcaller\\tstatus\\tn_bams\\n${runID}\\tDECoN\\t%s\\t%s\\n" "\$STATUS" "\$N_BAMS" \\
        > ${runID}_decon_caller_status.tsv
    [ "\$STATUS" = "ok" ] || echo "WARNING: DECoN status for ${runID}: \$STATUS (see ${runID}_decon.log)" >&2
    """
}

process RUN_CLEARCNV {

    tag "$runID"

    publishDir "${params.outdir}/runs/${runID}/02_clearcnv", mode: 'copy'

    input:
    tuple val(runID), path(retained_bams)

    output:
    tuple val(runID), path("${runID}_clearcnv_cnv_calls_final.tsv"), emit: calls
    path "${runID}_clearcnv_caller_status.tsv",                     emit: status
    path "${runID}_clearcnv.log"

    script:
    """
    CONFIG=${params.clearCNV_params_nf}
    OUTDIR="${runID}"
    mkdir -p "\$OUTDIR"
    STATUS="ok"
    N_BAMS=\$(grep -c . ${retained_bams} || true)

    export PATH="\$HOME/.conda/envs/mamba-env/bin:\$PATH"

    if ! clearCNV workflow_cnv_calling \\
        -w "\$OUTDIR" \\
        -p "${runID}" \\
        -r \$(yq -r '.fasta_file' \$CONFIG) \\
        -b ${retained_bams} \\
        -d ${params.caller_bed} \\
        -k \$(yq -r '.blacklist' \$CONFIG) \\
        -c \$(yq -r '.cores' \$CONFIG) \\
        --expected_artefacts \$(yq -r '.expected_artefacts' \$CONFIG) \\
        --sample_score_factor \$(yq -r '.sample_score_factor' \$CONFIG) \\
        --minimum_group_sizes \$(yq -r '.minimum_group_sizes' \$CONFIG) \\
        --zscale \$(yq -r '.zscale' \$CONFIG) \\
        --size \$(yq -r '.size' \$CONFIG) \\
        --del_cutoff \$(yq -r '.del_cutoff' \$CONFIG) \\
        --dup_cutoff \$(yq -r '.dup_cutoff' \$CONFIG) \\
        --trans_prob \$(yq -r '.trans_prob' \$CONFIG) > ${runID}_clearcnv.log 2>&1; then
        STATUS="failed"
    fi

    RESULTS="\$OUTDIR/${runID}/results/cnv_calls.tsv"
    if [ -f "\$RESULTS" ]; then
        cp "\$RESULTS" ${runID}_clearcnv_cnv_calls_final.tsv
    else
        [ "\$STATUS" = "ok" ] && STATUS="no_calls_file"
        : > ${runID}_clearcnv_cnv_calls_final.tsv
    fi

    printf "run\\tcaller\\tstatus\\tn_bams\\n${runID}\\tclearCNV\\t%s\\t%s\\n" "\$STATUS" "\$N_BAMS" \\
        > ${runID}_clearcnv_caller_status.tsv
    [ "\$STATUS" = "ok" ] || echo "WARNING: clearCNV status for ${runID}: \$STATUS (see ${runID}_clearcnv.log)" >&2
    """
}

// ----------------------------------------------------------------------
// Stages 3-5: framework filters
// ----------------------------------------------------------------------

process PRELIM_FILTER {

    tag "$runID"
    label 'r_small'

    publishDir "${params.outdir}/runs/${runID}/03_prelim_filter", mode: 'copy'

    input:
    tuple val(runID), path(decon_calls), path(clearcnv_calls), path(counts), path(outliers)

    output:
    tuple val(runID), path("${runID}_stage2_calls.tsv")

    script:
    """
    module load R
    prelim_filter.R --run ${runID} \\
        --decon ${decon_calls} --clearcnv ${clearcnv_calls} \\
        --counts ${counts} --outliers ${outliers} --bed ${params.framework_bed} \\
        --bf_min ${params.bf_min} --min_rc ${params.min_read_count} \\
        --utils ${utils}
    """
}

process CALLER_CONCORDANCE {

    tag "$runID"
    label 'r_small'

    publishDir "${params.outdir}/runs/${runID}/04_prioritisation", mode: 'copy'

    input:
    tuple val(runID), path(stage2)

    output:
    tuple val(runID), path("${runID}_stage3_calls.tsv")

    script:
    """
    module load R
    caller_concordance.R --run ${runID} --stage2 ${stage2} \\
        --reciprocal ${params.reciprocal_overlap} --utils ${utils}
    """
}

process ASSESS_RR_VAF {

    tag "$runID"
    label 'r_small'

    publishDir "${params.outdir}/runs/${runID}/05_rr_vaf", mode: 'copy'

    input:
    tuple val(runID), path(stage3), path(counts), path(outliers), val(vcf_dir)

    output:
    tuple val(runID), path("${runID}_stage4_calls.tsv"), emit: stage4
    path "${runID}_rr_per_exon.tsv",                    emit: rr
    path "${runID}_vaf_per_snv.tsv",                    emit: vaf

    script:
    """
    module load R
    assess_rr_vaf.R --run ${runID} --stage3 ${stage3} \\
        --counts ${counts} --outliers ${outliers} --vcf_dir "${vcf_dir}" \\
        --vcf_pattern '${params.vcf_pattern}' --id_regex '${params.sample_id_regex}' \\
        --n_ref ${params.n_ref} --rr_tol ${params.rr_tol} \\
        --rr_del_max ${params.rr_del_max} --rr_dup_min ${params.rr_dup_min} \\
        --rr_max_spread ${params.rr_max_spread} \\
        --vaf_max ${params.vaf_max} --dp_min ${params.dp_min} \\
        --dup_min_prop ${params.dup_min_prop_concordant} --dup_max_disc ${params.dup_max_discordant} \\
        --utils ${utils}
    """
}

// ----------------------------------------------------------------------
// Plots
// ----------------------------------------------------------------------

process PLOT_CALLS {

    tag "$runID"
    label 'r_small'

    publishDir "${params.outdir}/plots", mode: 'copy'

    input:
    tuple val(runID), path(stage4), path(stage2), path(counts), path(outliers), val(vcf_dir)

    output:
    path "calls/**", optional: true

    script:
    def exon_rds = params.exon_track_rds ?: ''
    """
    module load R
    plot_cnv_calls.R --run ${runID} --stage4 ${stage4} --stage2 ${stage2} \\
        --counts ${counts} --outliers ${outliers} --vcf_dir "${vcf_dir}" \\
        --bed ${params.framework_bed} --exon_track_rds '${exon_rds}' \\
        --vcf_pattern '${params.vcf_pattern}' --id_regex '${params.sample_id_regex}' \\
        --n_ref ${params.n_ref} --vaf_max ${params.vaf_max} --dp_min ${params.dp_min} \\
        --rr_min_max ${params.plot_rr_min_max} \\
        --show_vaf ${params.plot_vaf} --show_vaf_legend ${params.plot_vaf_legend} \\
        --show_raw_depth ${params.plot_raw_depth} \\
        --plot_prelim_rejected ${params.plot_prelim_rejected} \\
        --utils ${utils}
    """
}

process PLOT_MDS_ALL_RUNS {

    label 'r_small'

    publishDir "${params.outdir}/plots/outliers", mode: 'copy'

    input:
    path mds_files

    output:
    path "MDS_all_runs*", optional: true

    script:
    """
    module load R pandoc
    plot_mds_all_runs.R --utils ${utils}
    """
}

// ----------------------------------------------------------------------
// Cross-run tracking
// ----------------------------------------------------------------------

process COLLATE_FRAMEWORK {

    label 'r_small'

    publishDir "${params.outdir}/framework", mode: 'copy'

    input:
    path tables

    output:
    path "0*/*"
    path "*.tsv"

    script:
    """
    module load R
    collate_framework.R --utils ${utils}
    """
}

// ----------------------------------------------------------------------
// Workflow
// ----------------------------------------------------------------------

workflow {

    runs_ch = Channel
        .fromPath(params.samplesheet, checkIfExists: true)
        .splitCsv(header: true)
        .map { row -> tuple(row.run_id.trim(), row.bam_dir.trim(), row.vcf_dir.trim()) }

    vcf_ch = runs_ch.map { runID, bam_dir, vcf_dir -> tuple(runID, vcf_dir) }

    // Stage 1
    counts_ch = COUNT_READS(runs_ch)
    DETECT_OUTLIERS(counts_ch)

    // Stage 2
    RUN_DECON(DETECT_OUTLIERS.out.retained.join(counts_ch))
    RUN_CLEARCNV(DETECT_OUTLIERS.out.retained)

    // Stage 3
    prelim_in = RUN_DECON.out.calls
        .join(RUN_CLEARCNV.out.calls)
        .join(counts_ch)
        .join(DETECT_OUTLIERS.out.outliers)
    stage2_ch = PRELIM_FILTER(prelim_in)

    // Stage 4
    stage3_ch = CALLER_CONCORDANCE(stage2_ch)

    // Stage 5
    ASSESS_RR_VAF(
        stage3_ch
            .join(counts_ch)
            .join(DETECT_OUTLIERS.out.outliers)
            .join(vcf_ch)
    )

    // Plots (on by default; --skip_plots to turn off)
    if (!params.skip_plots) {
        PLOT_CALLS(
            ASSESS_RR_VAF.out.stage4
                .join(stage2_ch)
                .join(counts_ch)
                .join(DETECT_OUTLIERS.out.outliers)
                .join(vcf_ch)
        )
        PLOT_MDS_ALL_RUNS(DETECT_OUTLIERS.out.mds.collect())
    }

    // Tracking across all runs
    tables_ch = DETECT_OUTLIERS.out.outliers.map { it[1] }
        .mix(stage2_ch.map { it[1] })
        .mix(ASSESS_RR_VAF.out.stage4.map { it[1] })
        .mix(ASSESS_RR_VAF.out.rr)
        .mix(ASSESS_RR_VAF.out.vaf)
        .mix(RUN_DECON.out.status)
        .mix(RUN_CLEARCNV.out.status)
        .collect()
    COLLATE_FRAMEWORK(tables_ch)
}
