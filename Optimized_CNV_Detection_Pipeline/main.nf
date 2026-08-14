#!/usr/bin/env nextflow

// ======================================================================
// CNV calling pipeline with ClearCNV and DECoN
// ======================================================================

params.base_out = "/fs03/vh83/projects/temp_dnascreen_copy/dnascreen/Optimized_CNV_Detection_Pipeline/output"

// ----------------------------------------------------------------------
// Processes
// ----------------------------------------------------------------------

// THE MERGE RESULTS STILL NOT WORKING!
// a couple things that might not work with merge_results. the cnv_calls files are diff between clearCNV and decon.
// i think i might just have a tool-specific table and collate them tgt in R

process run_clearcnv {

    tag "$runID"

    input:
    tuple val(runID), path(bam_dir)

    output:
    tuple val(runID), path("${runID}_clearcnv_cnv_calls_final.tsv"), val("clearcnv")

    publishDir "${params.base_out}/clearcnv_output/${runID}", mode: 'copy'

    script:
    """
    OUTDIR="$runID"
    mkdir -p "\$OUTDIR/bams"

    BAM_LIST="\$OUTDIR/bams/bams_${runID}.txt"
    ls "$bam_dir"/*.bam | xargs -I{} realpath {} > "\$BAM_LIST"

    source activate mamba-env

    clearCNV workflow_cnv_calling \\
        -w "\$OUTDIR" \\
        -p "$runID" \\
        -r \$(yq -r '.fasta_file' ${params.clearCNV_params_nf}) \\
        -b "\$BAM_LIST" \\
        -d \$(yq -r '.bed_file' ${params.clearCNV_params_nf}) \\
        -k \$(yq -r '.blacklist' ${params.clearCNV_params_nf}) \\
        -c \$(yq -r '.cores' ${params.clearCNV_params_nf}) \\
        --expected_artefacts \$(yq -r '.expected_artefacts' ${params.clearCNV_params_nf}) \\
        --sample_score_factor \$(yq -r '.sample_score_factor' ${params.clearCNV_params_nf}) \\
        --minimum_group_sizes \$(yq -r '.minimum_group_sizes' ${params.clearCNV_params_nf}) \\
        --zscale \$(yq -r '.zscale' ${params.clearCNV_params_nf}) \\
        --size \$(yq -r '.size' ${params.clearCNV_params_nf}) \\
        --del_cutoff \$(yq -r '.del_cutoff' ${params.clearCNV_params_nf}) \\
        --dup_cutoff \$(yq -r '.dup_cutoff' ${params.clearCNV_params_nf}) \\
        --trans_prob \$(yq -r '.trans_prob' ${params.clearCNV_params_nf})

    cp "\$OUTDIR/$runID/results/cnv_calls.tsv" "${runID}_clearcnv_cnv_calls_final.tsv"
    """
}

process run_decon {

    tag "$runID"

    errorStrategy 'ignore'   // 👈 ignore errors and continue with other runs

    input:
    tuple val(runID), path(bam_dir)

    output:
    tuple val(runID), path("${runID}_decon_cnv_calls_final.tsv"), val("decon")

    publishDir "${params.base_out}/decon_output/${runID}", mode: 'copy'

    script:
    """
    module load R

    CONFIG=${params.decon_params_nf}
    RUNID="${runID}"   # make runID available in bash

    # --- Load YAML values ---
    BED=\$(yq -r '.bed_file' \$CONFIG)
    FASTA=\$(yq -r '.fasta_file' \$CONFIG)
    DECON_FOLDER=\$(yq -r '.deconFolder' \$CONFIG)
    MINCORR=\$(yq -r '.mincorr' \$CONFIG)
    MINCOV=\$(yq -r '.mincov' \$CONFIG)
    TRANSPROB=\$(yq -r '.transProb' \$CONFIG)

    # --- Prepare output dirs ---
    OUTDIR="\$RUNID"
    mkdir -p "\$OUTDIR/bams"

    BAM_TXT="\$OUTDIR/bams/bams_\${RUNID}.txt"
    ls "$bam_dir"/*.bam > "\$BAM_TXT"

    if [ ! -s "\$BAM_TXT" ]; then
        echo "⚠️ No BAM files found in $bam_dir"
        exit 1
    fi

    echo "Starting DECoN for runID \$RUNID at \$(date)"

    OUTPUT_BAMS="\$OUTDIR/output.bams"
    OUTPUT_RDATA="\$OUTDIR/output.bams.RData"

    echo "Running ReadInBams.R..."
    Rscript "\$DECON_FOLDER/ReadInBams.R" --bams "$bam_dir" --bed "\$BED" --fasta "\$FASTA" --out "\$OUTPUT_BAMS"

    echo "Running IdentifyFailures.R..."
    Rscript "\$DECON_FOLDER/IdentifyFailures.R" --RData "\$OUTPUT_RDATA" --mincorr "\$MINCORR" --mincov "\$MINCOV" --out "\$OUTDIR/failures"

    echo "Running makeCNVcalls.R..."
    Rscript "\$DECON_FOLDER/makeCNVcalls.R" --RData "\$OUTPUT_RDATA" --transProb "\$TRANSPROB" --plot None --out "\$OUTDIR" --failures "\$OUTDIR/failures_Failures.txt"

    # --- Find calls_all.txt anywhere inside OUTDIR ---
    CNV_FILE=\$(find . -name "*_all.txt" -type f | head -n 1)

    if [ -z "\$CNV_FILE" ]; then
        echo "⚠️ CNV calls file not found in current directory"
        exit 1
    fi

    cp "\$CNV_FILE" "${runID}_decon_cnv_calls_final.tsv"
    """
}

// Not using GATK anymore!
// process run_gatk {

//     tag "$runID"

//     errorStrategy 'ignore'   // 👈 ignore errors and continue with other runs

//     input:
//     tuple val(runID), path(bam_dir)

//     output:
//     tuple val(runID), path("${runID}_gatk_all_cnv_calls.txt"), val("gatk")

//     publishDir "${params.base_out}/gatk_output/${runID}", mode: 'copy'

//     script:
//     """
//     # Prepare params.yaml dynamically
//     CONFIG="${params.gatk_params_nf}"
//     OUTDIR="${runID}"

//     # Make OUTDIR absolute
//     mkdir -p "\$OUTDIR"
//     OUTDIR_ABS=\$(realpath "\$OUTDIR")
//     BAM_DIR_ABS=\$(realpath "${bam_dir}")

//     # Copy config template into run-specific folder
//     cp "\$CONFIG" "\$OUTDIR_ABS/params.yaml"

//     # Patch params.yaml with run-specific absolute BAM dir + absolute output folder + project name
//     yq -i '.bams_dir = "'"\$BAM_DIR_ABS"'"' "\$OUTDIR_ABS/params.yaml"
//     yq -i '.outputFolder = "'"\$OUTDIR_ABS"'"' "\$OUTDIR_ABS/params.yaml"
//     yq -i '.project_name = "'"${runID}"'"' "\$OUTDIR_ABS/params.yaml"

//     # Run the GATK CNV pipeline
//     bash "${params.gatk_main}" "\$OUTDIR_ABS/params.yaml"

//     # Final output = all_cnv_calls.txt renamed with runID
//     cp "\$OUTDIR_ABS/all_cnv_calls.txt" "${runID}_gatk_all_cnv_calls.txt"
//     """
// }




process merge_results {

    tag "merge_${tool}"

    input:
    tuple val(runID), path(run_output), val(tool)

    output:
    path "output/combined_${tool}_calls.tsv"

    script:
    """
    mkdir -p output

    # Initialize output file if not exists
    if [ ! -f output/combined_${tool}_calls.tsv ]; then
        # add RunID column at the front
        awk 'BEGIN{FS=OFS="\\t"} NR==1{print "RunID",\$0; next} {print runID,\$0}' runID=$runID "$run_output" \
            > output/combined_${tool}_calls.tsv
    else
        # skip header for subsequent files
        awk 'BEGIN{FS=OFS="\\t"} NR>1{print runID,\$0}' runID=$runID "$run_output" \
            >> output/combined_${tool}_calls.tsv
    fi
    """
}

// ----------------------------------------------------------------------
// Workflow definition
// ----------------------------------------------------------------------

workflow {

    bam_dirs_ch = Channel
        .fromPath(params.bam_lst)
        .flatMap { file -> file.text.readLines() }
        .map { it.trim() }
        .filter { it }
        .map { line ->
            def dir = file(line)
            tuple(dir.parent.baseName, dir)
        }

    results_ch = Channel.empty()

    if (params.run_clearcnv) {
        cnv_results_ch = run_clearcnv(bam_dirs_ch)
        results_ch = results_ch.mix(cnv_results_ch)
    }

    if (params.run_decon) {
        decon_results_ch = run_decon(bam_dirs_ch)
        results_ch = results_ch.mix(decon_results_ch)
    }

    if (params.run_gatk) {
        gatk_results_ch = run_gatk(bam_dirs_ch)
        results_ch = results_ch.mix(gatk_results_ch)
    }

    // results_ch.collect() | merge_results
}
