# DNAscreen_CNV

Scripts for the PhD thesis:

> **Detection of Germline Copy Number Variation in a Small Targeted-Sequencing Panel: Application for Population Genomic Screening**
> Zhu Jun (Joshua) Law, School of Translational Medicine, Monash University (2026)

DNA Screen was a pilot population genomic screening study. It sequenced 10,263 participants with a nine-gene targeted short-read panel covering Hereditary Breast and Ovarian Cancer (*BRCA1*, *BRCA2*, *PALB2*), Lynch syndrome (*MLH1*, *MSH2*, *MSH6*) and Familial Hypercholesterolaemia (*LDLR*, *APOB*, *PCSK9*). The original study reported SNVs and indels only. This project tested whether germline copy number variants (CNVs) can be detected reliably in the 9,012 BAMs from participants who consented to further research.

## Thesis aims and repository layout

The thesis hypothesis was that clinically relevant CNVs can be detected from the DNA Screen panel using an optimised read depth-based framework, and that targeted long-read sequencing can resolve their breakpoints and orientation well enough to improve clinical interpretation. Each aim has its own folder, matching a results chapter of the thesis:

| Folder | Thesis aim | Methods | Results chapter |
|---|---|---|---|
| [`Aim_1_scripts/`](Aim_1_scripts/) | **Aim 1.** Evaluate VS-CNV for CNV detection in DNA Screen, and build an aCGH-validated truth set | §2.2 | Chapter 3: Evaluation of VS-CNV in the DNA Screen panel |
| [`Aim_2_scripts/`](Aim_2_scripts/) | **Aim 2.** Benchmark alternative read depth-based CNV callers, and build an optimised CNV detection framework | §2.3 | Chapter 4: An optimised CNV detection framework for DNA Screen |
| [`Aim_3_scripts/`](Aim_3_scripts/) | **Aim 3.** Use targeted long-read sequencing (ONT and PacBio) to resolve breakpoints and duplication orientation | §2.4 | Chapter 5: Targeted long-read sequencing of DNA Screen CNVs |

The tables below list scripts in the order of the thesis Methods (Chapter 2), with the section that cites each one.

---

## Aim 1: Evaluation of VS-CNV (Chapter 3)

VS-CNV (VarSeq v2.6.0) was run per sequencing run (32 runs) at "Very High Sensitivity". In-silico single-exon CNVs and ClinVar pathogenic/likely pathogenic (PLP) CNVs were simulated to measure recall and to derive filtering thresholds. Prioritised calls were then validated by aCGH. Result: 10 confirmed CNVs, giving a truth set of 10 CNV-positive and 64 likely copy-neutral BAMs.

| Section | Step | Script(s) |
|---|---|---|
| 2.2.2 | Collate VS-CNV CNV- and BAM-level exports into a per-BAM, per-CNV table | [`vscnv_output_processing/summarise_CNVs.ipynb`](Aim_1_scripts/vscnv_output_processing/summarise_CNVs.ipynb) |
| 2.2.2 | Raw p < 0.05 filter, Benjamini–Hochberg adjustment, metric distributions | [`plots_on_vscnv_output/master_processing_real_data.R`](Aim_1_scripts/plots_on_vscnv_output/master_processing_real_data.R) |
| 2.2.2 | Z-score vs read-ratio scatter plot of all retained calls | [`plots_on_vscnv_output/z_score_rr_scatter_plot.R`](Aim_1_scripts/plots_on_vscnv_output/z_score_rr_scatter_plot.R) |
| 2.2.3 | Distribution of CNV calls per BAM | [`outlier_detection/summary_of_samples.R`](Aim_1_scripts/outlier_detection/summary_of_samples.R) |
| 2.2.3 | RPKM-normalised exon read depth → Pearson correlation → MDS | [`outlier_detection/correlation_matrix.R`](Aim_1_scripts/outlier_detection/correlation_matrix.R) |
| 2.2.3 | DBSCAN outlier BAM detection (minPts = 5, run-specific ε from the 5-NN elbow) | [`outlier_detection/automated_outlier_detection.R`](Aim_1_scripts/outlier_detection/automated_outlier_detection.R) |
| 2.2.4 | Exon read depths via ExomeDepth `getBamCounts`, cached as RData | [`PCA_ReadDepth/iterate_RD_all_runs.R`](Aim_1_scripts/PCA_ReadDepth/iterate_RD_all_runs.R) |
| 2.2.4 | PCA of normalised read depth, coloured by run | [`PCA_ReadDepth/PCA_RD_BAMS_diff_runs_gradient.R`](Aim_1_scripts/PCA_ReadDepth/PCA_RD_BAMS_diff_runs_gradient.R) |
| 2.2.5 | K-means on mean PC1/PC2 per run to choose cluster-representative runs for simulation | [`PCA_ReadDepth/mean_pca_find_distinct_runs.R`](Aim_1_scripts/PCA_ReadDepth/mean_pca_find_distinct_runs.R) |
| 2.2.6 | Single-exon CNV simulation (CN = 0–4, 165 exons) by read up/down-sampling; re-simulation of missed events | [`exon_sims/singleexon_sims/master_script/`](Aim_1_scripts/exon_sims/singleexon_sims/master_script/): `simulation_singlereg_master.sh`, `simulation_singlereg_missed_master.sh`, `simulation_repeat_cn{0..4}.sh` |
| 2.2.7 | ClinVar PLP CNV simulation (533 deletions at CN = 1, 88 duplications at CN = 3, in 1,582 BAMs from runs 6, 13, 14, 18, 29, 33) | [`exon_sims/clinvar_sims/master_script/`](Aim_1_scripts/exon_sims/clinvar_sims/master_script/): `simulation_multiregion_master_clinvar.sh`, `simulation_multiregion_missed_master_clinvar.sh` |
| 2.2.8 | Match VS-CNV calls on simulated BAMs to simulated regions; compute recall | [`vscnv_output_processing/collating_sim_CNV_singlereg.ipynb`](Aim_1_scripts/vscnv_output_processing/collating_sim_CNV_singlereg.ipynb), [`collating_sim_CNV_clinvar.ipynb`](Aim_1_scripts/vscnv_output_processing/collating_sim_CNV_clinvar.ipynb) |
| 2.2.8 | Simulation summaries and plots | [`simulation_cnv_plots/master_processing_singleexon_sims.R`](Aim_1_scripts/simulation_cnv_plots/master_processing_singleexon_sims.R), [`master_processing_clinvar_sims.R`](Aim_1_scripts/simulation_cnv_plots/master_processing_clinvar_sims.R) |
| 2.2.8 | Percentage of exonic overlap for missed ClinVar CNVs | [`simulation_cnv_plots/clinvar_cnvs_exon_percentage_overlap_analysis.R`](Aim_1_scripts/simulation_cnv_plots/clinvar_cnvs_exon_percentage_overlap_analysis.R) |
| 2.2.9 | Apply simulation-derived z-score and read-ratio thresholds to real calls | [`plots_on_vscnv_output/clinvar_naive_cnvs_volcano_plots.R`](Aim_1_scripts/plots_on_vscnv_output/clinvar_naive_cnvs_volcano_plots.R) |
| 2.2.9, 2.2.10 | Selection of calls for aCGH; z-score and read ratio of aCGH-confirmed vs unsupported calls (AUROC, Wilson CIs) | [`plots_on_vscnv_output/compare_metrics_with_acgh_confirmation.R`](Aim_1_scripts/plots_on_vscnv_output/compare_metrics_with_acgh_confirmation.R) |
| — | Shared helpers and supplementary plots | [`plots_on_vscnv_output/functions_for_processing_sim_real_calls.R`](Aim_1_scripts/plots_on_vscnv_output/functions_for_processing_sim_real_calls.R), [`summary_of_cnv_calls.R`](Aim_1_scripts/plots_on_vscnv_output/summary_of_cnv_calls.R), [`per_exon_RD_distribution.R`](Aim_1_scripts/plots_on_vscnv_output/per_exon_RD_distribution.R) |

The ClinVar PLP input files named in §2.2.7 (`clinvar_plp_cn1.bed`, `clinvar_plp_cn3.bed`) are not included in this repository. Scripts under `deprecated_scripts/` or with `_deprecated`/`_old` suffixes are earlier versions kept for reference; they were not used for the thesis results.

---

## Aim 2: Optimised CNV detection framework (Chapter 4)

Eleven germline CNV callers plus VS-CNV were benchmarked on the truth set. DECoN gave the highest recall and ClearCNV the highest precision. Both were built into a per-run framework that combines caller concordance with checks on exon-level read depth, read ratio and VAF. Applied to the full dataset, the framework found 12 high-confidence CNVs, which were then clinically curated to ACMG/ClinGen standards.

### [`CNVbenchmarkeR2/`](Aim_2_scripts/CNVbenchmarkeR2/): caller benchmarking (§2.3.1–2.3.4)

This is a modified copy of [CNVbenchmarkeR2](https://github.com/jpuntomarcos/CNVbenchmarkeR2) (Munté et al., 2024). Per-tool wrappers and parameters are in [`tools/`](Aim_2_scripts/CNVbenchmarkeR2/tools/): ClearCNV, GATK-gCNV (`germlineCNVcaller`), Atlas-CNV, Cobalt, ClinCNV, CNVkit, VisCap, DECoN, panelcn.MOPS, ExomeDepth, CODEX2 (CoNVaDING is included but could not be run on the HPC). Run it with `runBenchmark.sh` / `runBenchmark.R`, configured through `datasets.yaml` and `tools.yaml`.

| Section | Step | Script(s) in `Rscripts/` |
|---|---|---|
| 2.3.1 | Classify calls as true/false positive (≥ 20% exon overlap, ≥ 50% of called exons concordant with aCGH); recall/precision with Wilson CIs | [`assign_acgh_concordance_and_check_nonvs_cnv.R`](Aim_2_scripts/CNVbenchmarkeR2/Rscripts/assign_acgh_concordance_and_check_nonvs_cnv.R) |
| 2.3.1 | Jitter plot of z-scaled per-tool quality metrics | [`jitter_plot.R`](Aim_2_scripts/CNVbenchmarkeR2/Rscripts/jitter_plot.R) |
| 2.3.2 | Shared and unique per-run DECoN/ClearCNV calls (Venn) | [`per_run_decon_clearcnv_venn.R`](Aim_2_scripts/CNVbenchmarkeR2/Rscripts/per_run_decon_clearcnv_venn.R) |
| 2.3.3–2.3.4 | Read depth, read ratio and VAF of true vs false positive DECoN calls; manual-inspection decision tree | [`optimized_framework_calls.R`](Aim_2_scripts/CNVbenchmarkeR2/Rscripts/optimized_framework_calls.R) |
| 2.3.3 | Exon-level read-ratio + VAF (Wilson CI) plots per call | [`plot_rr_with_vaf.R`](Aim_2_scripts/CNVbenchmarkeR2/Rscripts/plot_rr_with_vaf.R) |
| — | Shared helpers | [`util_functions_for_tool_cnvcalls.R`](Aim_2_scripts/CNVbenchmarkeR2/Rscripts/util_functions_for_tool_cnvcalls.R) |

The other scripts in `Rscripts/` (`*_hgsa.R`, `*_run31.R`, `*_archived.R`, the violin and accuracy plots) are exploratory or conference analyses and were not cited in the thesis.

### [`Optimized_CNV_Detection_Pipeline/`](Aim_2_scripts/Optimized_CNV_Detection_Pipeline/): the framework (§2.3.2, §2.3.5)

This is a Nextflow DSL2 pipeline ([`main.nf`](Aim_2_scripts/Optimized_CNV_Detection_Pipeline/main.nf), [`nextflow.config`](Aim_2_scripts/Optimized_CNV_Detection_Pipeline/nextflow.config)) that runs the framework described in §2.3.5 and Figure 4.7 on each sequencing run:

1. Outlier BAM removal (ExomeDepth counts → RPKM correlation → MDS → DBSCAN; §2.2.3)
2. Per-run CNV calling with DECoN and ClearCNV ([`decon/decon_params_nf.yaml`](Aim_2_scripts/Optimized_CNV_Detection_Pipeline/decon/decon_params_nf.yaml), [`clearCNV/clearCNV_params_nf.yaml`](Aim_2_scripts/Optimized_CNV_Detection_Pipeline/clearCNV/clearCNV_params_nf.yaml))
3. Preliminary filter: DECoN Bayes factor ≥ 10; mean read depth ≥ 100 for single-exon calls
4. Caller concordance: calls sharing ≥ 50% of exons are high confidence; the rest go to manual inspection
5. Read-ratio and VAF decision tree for single-caller calls (§2.3.4, Figure 2.5)

Usage, samplesheet format and outputs are documented in the [pipeline README](Aim_2_scripts/Optimized_CNV_Detection_Pipeline/README.md).

---

## Aim 3: Targeted long-read sequencing (Chapter 5)

Samples carrying the 12 high-confidence CNVs were enriched with a custom hybridisation panel and sequenced on ONT PromethION 2 (R10.4.1) and PacBio Vega. Breakpoints were resolved to single-base-pair resolution by manually inspecting split reads in IGV. This confirmed two VUS duplications as tandem, so they were reclassified as pathogenic.

| Section | Step | Script(s) |
|---|---|---|
| 2.4.2 | ONT: Dorado basecalling and demultiplexing, read-group tagging, duplicate marking, per-barcode HsMetrics, MultiQC (two sequencing rounds) | [`ont_demultiplexing_round1/`](Aim_3_scripts/ont_demultiplexing_round1/), [`ont_demultiplexing_round2/`](Aim_3_scripts/ont_demultiplexing_round2/) |
| 2.4.2 | PacBio: HiFiTargetEnrichment Snakemake workflow, modified to realign with minimap2 instead of pbmm2 | [`pacbio_demultiplexing/workflow_HiFiTargetEnrichment/`](Aim_3_scripts/pacbio_demultiplexing/workflow_HiFiTargetEnrichment/) (see its [README](Aim_3_scripts/pacbio_demultiplexing/workflow_HiFiTargetEnrichment/README.MD)) |
| 2.4.2 | QC metrics: read length, Phred quality, % target bases ≥ 100X | `multiqc.sh` (in each processing folder above) |
| 2.4.2 | SV calling with Sniffles2, CuteSV, SVIM, DeBreak (default parameters) | [`ont_pacbio_sv_callers/{sniffles,cutesv,svim,debreak}/*_loop.sh`](Aim_3_scripts/ont_pacbio_sv_callers/) (`*_round2.sh` for the second ONT round) |
| 2.4.2 | Concordance of SV calls with aCGH-confirmed CNVs (same type, same chromosome, ≥ 1 bp overlap) | [`ont_pacbio_sv_callers/master_table.py`](Aim_3_scripts/ont_pacbio_sv_callers/master_table.py), [`SV_tool_recall_with_acgh.R`](Aim_3_scripts/SV_tool_recall_with_acgh.R) |
| 2.4.3 | Count reads supporting each manually resolved junction sequence | [`check_supporting_reads.sh`](Aim_3_scripts/check_supporting_reads.sh) |
| 2.4.5 | Total absolute breakpoint offset of SV-caller calls vs manually resolved breakpoints | [`ont_pacbio_sv_callers/{sniffles,cutesv,svim,debreak}/compare_*_cnv.py`](Aim_3_scripts/ont_pacbio_sv_callers/) |
| 2.4.5 | Offset plot per BAM and platform | [`script_for_offset_plotting.R`](Aim_3_scripts/script_for_offset_plotting.R) |
| Ch. 5 | Short-read vs long-read CNV interval comparison with exon/transcript track | [`breakpoint_uncertainty.R`](Aim_3_scripts/breakpoint_uncertainty.R) |
| — | Helpers: VCF → BED filtering, per-call read support, call-threshold exploration | [`ont_pacbio_sv_callers/filter_vcfs_to_bed.py`](Aim_3_scripts/ont_pacbio_sv_callers/filter_vcfs_to_bed.py), [`read_support_for_vcf_calls.py`](Aim_3_scripts/ont_pacbio_sv_callers/read_support_for_vcf_calls.py), [`find_zero_call_thresholds.py`](Aim_3_scripts/ont_pacbio_sv_callers/find_zero_call_thresholds.py), [`find_shared_calls_bc1001_bc1099.py`](Aim_3_scripts/ont_pacbio_sv_callers/find_shared_calls_bc1001_bc1099.py) |

Manual breakpoint resolution in IGV (§2.4.3) and clinical re-curation (§2.4.4) were done by hand and have no scripts.

---

## Main software

| Tool | Version | Used in |
|---|---|---|
| VarSeq / VS-CNV | 2.6.0 | Aim 1, benchmarking |
| R | 4.2.3 | Aims 1–3 |
| Python | 3.9.23 | Aims 1, 3 |
| samtools / bedtools | 1.19.3 / 2.31.1 | CNV simulation, read support |
| ExomeDepth | 1.1.17 | Read counts, benchmarking |
| dbscan (R) | 1.2 | Outlier detection |
| DECoN / ClearCNV | 2.0.1 / 0.306 | Framework |
| VariantAnnotation (R) | 1.44.1 | VAF extraction |
| Nextflow | 24.04.3 | Framework pipeline |
| Dorado | 1.2.0 | ONT basecalling |
| minimap2, Sniffles2, CuteSV, SVIM, DeBreak | see thesis Fig. 2.6 | Aim 3 |
| MultiQC / Picard | 1.23 / 2.19.0 | Long-read QC |
| IGV | 2.19.8 | Manual breakpoint inspection |

Versions for the other benchmarked callers are listed in thesis Table 2.2. Most scripts were written for SLURM on the Monash M3 HPC cluster and contain hard-coded paths, which will need editing before use elsewhere.

## License

See [LICENSE](LICENSE). The bundled CNVbenchmarkeR2 and PacBio HiFiTargetEnrichment workflows keep their original licenses.
