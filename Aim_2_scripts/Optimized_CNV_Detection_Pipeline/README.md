# Optimized_CNV_Detection_Pipeline

Nextflow implementation of the optimised CNV framework, applied per sequencing run.

| Stage | Process | What it does |
|---|---|---|
| 1 | `COUNT_READS`, `DETECT_OUTLIERS` | ExomeDepth counts on the framework BED; RPKM → sample correlation → MDS → DBSCAN. DBSCAN noise and zero-coverage BAMs are removed |
| 2 | `RUN_DECON`, `RUN_CLEARCNV` | CNV calling on retained BAMs only |
| 3 | `PRELIM_FILTER` | DECoN BF ≥ 10 (all calls); mean read depth ≥ 100 (single-exon calls only, both callers) |
| 4 | `CALLER_CONCORDANCE` | Same sample + type, reciprocal overlap ≥ 0.5 → dual-caller (high confidence); otherwise single-caller |
| 5 | `ASSESS_RR_VAF` | Single-caller calls: RR1/RR2/RR3 then the VAF decision tree |
| — | `PLOT_CALLS`, `PLOT_MDS_ALL_RUNS` | Per-call RR + VAF plots, outlier MDS plots (on by default) |
| — | `COLLATE_FRAMEWORK` | Cross-run tracking tables |

## Run

```bash
cd Optimized_CNV_Detection_Pipeline
sbatch run_framework.slurm                       # runs 4-35 (samplesheet_runs4-35.csv)
sbatch run_framework.slurm --samplesheet one_run.csv --outdir output/test_run7
sbatch run_framework.slurm -resume
```

Samplesheet (CSV): `run_id,bam_dir,vcf_dir`. VCFs are matched to BAMs by file stem, then by `DNS\d+` ID.
All thresholds are params in `nextflow.config`. `--skip_plots` turns plotting off, and `--plot_prelim_rejected true` also plots calls that fail stage 3.
`caller_bed` overrides `bed_file` in the DECoN and clearCNV YAMLs.

## Outputs (`output/optimised_framework/`)

```
framework/
  01_outliers/        outlier_bams.tsv, retained_bams.tsv
  02_raw_calls/       decon_calls.tsv, clearcnv_calls.tsv
  03_prelim_filter/   passed.tsv, failed.tsv (prelim_fail_reason)
  04_prioritisation/  dual_caller_high_confidence.tsv, dual_caller_clearcnv_partners.tsv,
                      single_caller_manual_inspection.tsv
  05_rr_vaf/          passed.tsv, failed_rr.tsv, failed_vaf.tsv, rr_per_exon.tsv, vaf_per_snv.tsv
  06_final/           high_confidence_calls.tsv, rejected_calls.tsv
  call_ledger.tsv     every raw call: status at each stage, fail_reason, final_status
  sample_ledger.tsv   every BAM: outlier status, calls reaching each stage
  framework_summary.tsv / framework_summary_by_run.tsv   counts per stage and caller
  run_status.tsv      caller status per run (ok / failed / no_calls_file)
plots/
  outliers/           <run>_MDS_plot.{png,html}, <run>_correlation_histogram.png, MDS_all_runs.html
  calls/<final_status>/<run>_<ID>_<gene>_<type>_<caller>_chr<c>-<start>-<end>.jpeg
runs/<run>/           per-run intermediate files
pipeline_info/        report.html, trace.tsv
```

`final_status` is one of: `rejected_prelim_filter`, `high_confidence_dual_caller`,
`merged_into_dual_caller` (the clearCNV partner of a dual-caller call),
`high_confidence_rr_vaf`, `rejected_rr`, `rejected_vaf`.
