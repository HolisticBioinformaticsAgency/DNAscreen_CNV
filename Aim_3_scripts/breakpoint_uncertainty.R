# ==============================================================
# Short-read vs long-read CNV interval comparison, with exon-
# level transcript track and terminal-exon UTR representation.
#
# Long-read intervals are displayed above short-read intervals.
#
# Short-read boxes may be widened locally for label legibility,
# without rescaling genomic coordinates used for transcript, exon,
# UTR, long-read, uncertainty, or tick positioning.
#
# Dotted guide lines link every displayed short-read interval edge
# to its true genomic breakpoint coordinate at the top of the exon
# track.
#
# IMPORTANT:
# Bounded short-read uncertainty ticks are anchored to the actual
# genomic edge of the adjacent exon, rather than reconstructed from
# the call Start/End coordinate and intron length.
# ==============================================================

library(dplyr)
library(readr)
library(stringr)
library(ggplot2)
library(scales)
library(forcats)
library(purrr)
library(tidyr)

# --------------------------- Paths ----------------------------

high_confidence_file <- "/Users/zlaw0001/Library/CloudStorage/OneDrive-MonashUniversity/Joshua_PhD_Project/Aim_2/Optimized_CNV_Detection_Pipeline/rr_of_12_highconf_calls_with_vaf/high_confidence_calls.csv"

junction_file <- "/Users/zlaw0001/Library/CloudStorage/OneDrive-MonashUniversity/Joshua_PhD_Project/Aim_3/breakpoint_identification/ont_pacbio_sv_callers/sa_junction_read_counts.tsv"

exon_rds_file <- "/Users/zlaw0001/Library/CloudStorage/OneDrive-MonashUniversity/Joshua_PhD_Project/Aim_3/breakpoint_identification/ont_pacbio_sv_callers/cached_transcript_with_UTR_data.rds"

output_prefix <- "breakpoint_precision"

# ---------------- Tunable display parameters -------------------

# Larger values reduce the minimum imposed display width, resulting
# in less widening for narrow short-read CNV boxes.
label_chars_per_window <- 140

# Decrease towards 1.0 for a tighter text fit.
label_padding_multiplier <- 0.7

# Full-height coding-exon box.
exon_track_y_min <- 0.16
exon_track_y_max <- 0.44

# Thin UTR box, centred on the transcript backbone.
utr_track_y_min <- 0.23
utr_track_y_max <- 0.37

# This function creates a rounded DISPLAY LABEL only.
# All genomic coordinate calculations remain in exact base pairs.
format_genomic_length_label <- function(length_bp) {
  paste0(round(length_bp / 1000, 1), " kb")
}

clinically_relevant_transcripts <- c(
  BRCA1 = "NM_007294",
  BRCA2 = "NM_000059",
  PALB2 = "NM_024675",
  MSH2 = "NM_000251",
  MSH6 = "NM_000179",
  MLH1 = "NM_000249",
  APOB = "NM_000384",
  PCSK9 = "NM_174936",
  LDLR = "NM_000527"
)

canonical_chromosomes <- c(as.character(1:22), "X", "Y", "MT")

# ------------------------ Read inputs --------------------------

hc_calls <- read_csv(
  high_confidence_file,
  show_col_types = FALSE
) %>%
  mutate(
    sample_id_trimmed = as.character(sample_id_trimmed),
    Gene = as.character(Gene),
    CNV_Type = str_to_title(as.character(CNV_Type)),
    Chromosome = str_remove(as.character(Chromosome), "^chr"),
    Start = as.numeric(Start),
    End = as.numeric(End)
  )

junctions <- read_tsv(
  junction_file,
  show_col_types = FALSE
) %>%
  mutate(
    sample_id = as.character(sample_id),
    resolved_longread_junction = !is.na(junction_seq) &
      junction_seq != "NA"
  ) %>%
  filter(resolved_longread_junction) %>%
  distinct(sample_id, .keep_all = TRUE) %>%
  select(sample_id, resolved_longread_junction)

resolved_calls <- hc_calls %>%
  inner_join(
    junctions,
    by = c("sample_id_trimmed" = "sample_id")
  )

if (nrow(resolved_calls) == 0) {
  stop(
    "No resolved CNVs were found after joining high-confidence calls ",
    "to sa_junction_read_counts.tsv."
  )
}

longread_breakpoints <- tribble(
  ~sample_id_trimmed, ~lrs_chrom, ~lrs_start, ~lrs_end,
  "DNS004937", "17", 43045048, 43046166,
  "DNS004024", "17", 43111616, 43117534,
  "DNS006916", "17", 43078282, 43084362,
  "DNS006946", "17", 43078282, 43084362,
  "DNS000601", "13", 32325651, 32327136,
  "DNS006443", "16", 23605575, 23615114,
  "DNS010349", "2", 20993431, 21026997,
  "DNS008387", "19", 11094251, 11175478,
  "DNS006683", "1", 55029215, 55043368
) %>%
  mutate(
    lrs_chrom = str_remove(as.character(lrs_chrom), "^chr"),
    lrs_size_bp = lrs_end - lrs_start + 1L
  )

# ---------------- Load transcript/exon/UTR data -----------------

cached_annotation <- readRDS(exon_rds_file)

if (!is.list(cached_annotation)) {
  stop(
    "The RDS file must contain the cached transcript annotation list. ",
    "Expected components include $exons and $utr_segments."
  )
}

if (!all(c("exons", "utr_segments") %in% names(cached_annotation))) {
  stop(
    "The RDS object is missing required components: $exons and/or ",
    "$utr_segments."
  )
}

exonscombined <- cached_annotation$exons
utr_segments_raw <- cached_annotation$utr_segments

required_exon_columns <- c(
  "gene",
  "chromosome_name",
  "exon_chrom_start",
  "exon_chrom_end",
  "refseq_mrna",
  "strand",
  "rank"
)

missing_exon_columns <- setdiff(
  required_exon_columns,
  colnames(exonscombined)
)

if (length(missing_exon_columns) > 0) {
  stop(
    "The exon table in the RDS object is missing: ",
    paste(missing_exon_columns, collapse = ", ")
  )
}

required_utr_columns <- c(
  "gene",
  "chromosome_name",
  "utr_chrom_start",
  "utr_chrom_end",
  "refseq_mrna",
  "strand",
  "exon_rank",
  "utr_type"
)

missing_utr_columns <- setdiff(
  required_utr_columns,
  colnames(utr_segments_raw)
)

if (length(missing_utr_columns) > 0) {
  stop(
    "The UTR table in the RDS object is missing: ",
    paste(missing_utr_columns, collapse = ", ")
  )
}

exons <- exonscombined %>%
  mutate(
    gene = as.character(gene),
    chromosome_name = str_remove(as.character(chromosome_name), "^chr"),
    exon_chrom_start = as.numeric(exon_chrom_start),
    exon_chrom_end = as.numeric(exon_chrom_end),
    refseq_mrna = as.character(refseq_mrna),
    strand = as.numeric(strand),
    rank = as.integer(rank)
  ) %>%
  filter(
    gene %in% names(clinically_relevant_transcripts),
    refseq_mrna == clinically_relevant_transcripts[gene],
    chromosome_name %in% canonical_chromosomes
  ) %>%
  arrange(gene, chromosome_name, exon_chrom_start, exon_chrom_end) %>%
  group_by(gene, chromosome_name) %>%
  mutate(
    exon_number_genomic = row_number(),
    n_exons_gene = n()
  ) %>%
  ungroup() %>%
  rename(
    chromosomename = chromosome_name,
    exonchromstart = exon_chrom_start,
    exonchromend = exon_chrom_end,
    exon_rank = rank
  ) %>%
  mutate(
    exon_number_bio = if_else(
      strand == -1,
      n_exons_gene - exon_number_genomic + 1L,
      exon_number_genomic
    )
  )

utr_segments <- utr_segments_raw %>%
  mutate(
    gene = as.character(gene),
    chromosome_name = str_remove(as.character(chromosome_name), "^chr"),
    refseq_mrna = as.character(refseq_mrna),
    strand = as.numeric(strand),
    exon_rank = as.integer(exon_rank),
    utr_type = as.character(utr_type),
    utr_chrom_start = as.numeric(utr_chrom_start),
    utr_chrom_end = as.numeric(utr_chrom_end)
  ) %>%
  filter(
    gene %in% names(clinically_relevant_transcripts),
    refseq_mrna == clinically_relevant_transcripts[gene],
    chromosome_name %in% canonical_chromosomes,
    utr_type %in% c("5_prime_UTR", "3_prime_UTR")
  ) %>%
  rename(
    chromosomename = chromosome_name,
    utr_start = utr_chrom_start,
    utr_end = utr_chrom_end
  ) %>%
  arrange(gene, chromosomename, utr_start, utr_end)

if (nrow(exons) == 0) {
  stop("No exons remained after filtering to clinically relevant transcripts.")
}

if (nrow(utr_segments) == 0) {
  warning(
    "No UTR segments remained after filtering. Exons will still be plotted, ",
    "but no thinner UTR boxes will be shown."
  )
}

# ----------------- Infer short-read uncertainty -----------------

infer_breakpoint_uncertainty <- function(
    sample_id,
    gene,
    chromosome,
    cnv_start,
    cnv_end,
    cnv_type,
    exon_table
) {
  gene_exons <- exon_table %>%
    filter(
      gene == !!gene,
      chromosomename == str_remove(as.character(chromosome), "^chr")
    ) %>%
    arrange(exonchromstart)
  
  if (nrow(gene_exons) == 0) {
    warning("No exon coordinates found for ", sample_id, " (", gene, ").")
    return(NULL)
  }
  
  called_exons <- gene_exons %>%
    filter(
      exonchromstart <= cnv_end,
      exonchromend >= cnv_start
    )
  
  if (nrow(called_exons) == 0) {
    warning(
      "No exons overlap short-read CNV coordinates for ",
      sample_id,
      " (",
      gene,
      ")."
    )
    return(NULL)
  }
  
  first_called <- min(called_exons$exon_number_genomic)
  last_called <- max(called_exons$exon_number_genomic)
  total_exons <- max(gene_exons$n_exons_gene)
  
  first_exon <- gene_exons %>%
    filter(exon_number_genomic == first_called)
  
  last_exon <- gene_exons %>%
    filter(exon_number_genomic == last_called)
  
  # Left uncertainty: from the previous exon edge to the first called exon.
  if (first_called == 1) {
    left_uncertainty <- tibble(
      sample_id_trimmed = sample_id,
      breakpoint_side = "left",
      uncertainty_bp = NA_real_,
      boundary_coordinate = NA_real_,
      bounded = FALSE,
      uncertainty_label = "Unbounded"
    )
  } else {
    preceding_exon <- gene_exons %>%
      filter(exon_number_genomic == first_called - 1L)
    
    intron_start <- preceding_exon$exonchromend + 1L
    intron_end <- first_exon$exonchromstart - 1L
    intron_length <- intron_end - intron_start + 1L
    
    left_uncertainty <- tibble(
      sample_id_trimmed = sample_id,
      breakpoint_side = "left",
      uncertainty_bp = intron_length,
      boundary_coordinate = preceding_exon$exonchromend,
      bounded = TRUE,
      uncertainty_label = format_genomic_length_label(intron_length)
    )
  }
  
  # Right uncertainty: from the final called exon to the next exon edge.
  if (last_called == total_exons) {
    right_uncertainty <- tibble(
      sample_id_trimmed = sample_id,
      breakpoint_side = "right",
      uncertainty_bp = NA_real_,
      boundary_coordinate = NA_real_,
      bounded = FALSE,
      uncertainty_label = "Unbounded"
    )
  } else {
    following_exon <- gene_exons %>%
      filter(exon_number_genomic == last_called + 1L)
    
    intron_start <- last_exon$exonchromend + 1L
    intron_end <- following_exon$exonchromstart - 1L
    intron_length <- intron_end - intron_start + 1L
    
    right_uncertainty <- tibble(
      sample_id_trimmed = sample_id,
      breakpoint_side = "right",
      uncertainty_bp = intron_length,
      boundary_coordinate = following_exon$exonchromstart,
      bounded = TRUE,
      uncertainty_label = format_genomic_length_label(intron_length)
    )
  }
  
  bind_rows(left_uncertainty, right_uncertainty)
}

precision_data <- pmap_dfr(
  resolved_calls %>%
    select(
      sample_id_trimmed,
      Gene,
      Chromosome,
      Start,
      End,
      CNV_Type
    ),
  function(sample_id_trimmed, Gene, Chromosome, Start, End, CNV_Type) {
    infer_breakpoint_uncertainty(
      sample_id = sample_id_trimmed,
      gene = Gene,
      chromosome = Chromosome,
      cnv_start = Start,
      cnv_end = End,
      cnv_type = CNV_Type,
      exon_table = exons
    )
  }
)

uncertainty_wide <- precision_data %>%
  pivot_wider(
    id_cols = sample_id_trimmed,
    names_from = breakpoint_side,
    values_from = c(
      uncertainty_bp,
      boundary_coordinate,
      bounded,
      uncertainty_label
    ),
    names_glue = "{breakpoint_side}_{.value}"
  )

# ------------------- Exon-range lookup function -----------------

get_exon_range <- function(
    gene,
    chromosome,
    interval_start,
    interval_end,
    exon_table
) {
  gene_exons <- exon_table %>%
    filter(
      gene == !!gene,
      chromosomename == str_remove(as.character(chromosome), "^chr")
    ) %>%
    arrange(exonchromstart)
  
  if (nrow(gene_exons) == 0) {
    return(
      tibble(
        exon_range_label = NA_character_,
        matched_exon_start = NA_real_,
        matched_exon_end = NA_real_,
        n_matched_exons = 0L
      )
    )
  }
  
  overlapping_exons <- gene_exons %>%
    filter(
      exonchromstart <= interval_end,
      exonchromend >= interval_start
    )
  
  if (nrow(overlapping_exons) == 0) {
    return(
      tibble(
        exon_range_label = NA_character_,
        matched_exon_start = NA_real_,
        matched_exon_end = NA_real_,
        n_matched_exons = 0L
      )
    )
  }
  
  first_exon <- min(overlapping_exons$exon_number_bio)
  last_exon <- max(overlapping_exons$exon_number_bio)
  
  exon_range_label <- if (first_exon == last_exon) {
    paste0("exon ", first_exon)
  } else {
    paste0("exons ", first_exon, "-", last_exon)
  }
  
  tibble(
    exon_range_label = exon_range_label,
    matched_exon_start = min(overlapping_exons$exonchromstart),
    matched_exon_end = max(overlapping_exons$exonchromend),
    n_matched_exons = nrow(overlapping_exons)
  )
}

lrs_exon_ranges <- longread_breakpoints %>%
  inner_join(
    resolved_calls %>%
      select(sample_id_trimmed, Gene) %>%
      distinct(),
    by = "sample_id_trimmed"
  ) %>%
  mutate(
    exon_info = pmap(
      list(Gene, lrs_chrom, lrs_start, lrs_end),
      function(Gene, lrs_chrom, lrs_start, lrs_end) {
        get_exon_range(
          gene = Gene,
          chromosome = lrs_chrom,
          interval_start = lrs_start,
          interval_end = lrs_end,
          exon_table = exons
        )
      }
    )
  ) %>%
  unnest(exon_info) %>%
  select(sample_id_trimmed, exon_range_label)

sr_exon_match <- resolved_calls %>%
  select(sample_id_trimmed, Gene, Chromosome, Start, End) %>%
  mutate(
    exon_info = pmap(
      list(Gene, Chromosome, Start, End),
      function(Gene, Chromosome, Start, End) {
        get_exon_range(
          gene = Gene,
          chromosome = Chromosome,
          interval_start = Start,
          interval_end = End,
          exon_table = exons
        )
      }
    )
  ) %>%
  unnest(exon_info) %>%
  select(
    sample_id_trimmed,
    matched_exon_start,
    matched_exon_end,
    n_matched_exons
  )

# ---------- Build one row per CNV: short-read + long-read -------

sample_order <- c(
  "DNS004937",
  "DNS004024",
  "DNS006916",
  "DNS006946",
  "DNS000601",
  "DNS006443",
  "DNS010349",
  "DNS008387",
  "DNS006683"
)

cnv_data <- resolved_calls %>%
  select(sample_id_trimmed, Gene, Chromosome, Start, End, CNV_Type) %>%
  inner_join(
    longread_breakpoints,
    by = "sample_id_trimmed"
  ) %>%
  left_join(
    lrs_exon_ranges,
    by = "sample_id_trimmed"
  ) %>%
  left_join(
    uncertainty_wide,
    by = "sample_id_trimmed"
  ) %>%
  left_join(
    sr_exon_match,
    by = "sample_id_trimmed"
  ) %>%
  mutate(
    left_bounded = coalesce(left_bounded, FALSE),
    right_bounded = coalesce(right_bounded, FALSE),
    
    left_uncertainty_label = coalesce(
      left_uncertainty_label,
      "Unbounded"
    ),
    
    right_uncertainty_label = coalesce(
      right_uncertainty_label,
      "Unbounded"
    ),
    
    cnv_midpoint = (lrs_start + lrs_end) / 2,
    sr_size_bp = End - Start + 1L,
    
    sample_id_trimmed = factor(
      sample_id_trimmed,
      levels = sample_order
    ),
    
    panel_label = paste0(
      sample_id_trimmed,
      "\n",
      Gene,
      "\n",
      coalesce(exon_range_label, "exon range unknown"),
      "\n",
      CNV_Type
    ),
    
    # True genomic coordinates.
    sr_xmin_true = Start - cnv_midpoint,
    sr_xmax_true = End - cnv_midpoint,
    
    lr_xmin_true = lrs_start - cnv_midpoint,
    lr_xmax_true = lrs_end - cnv_midpoint,
    
    # Labels inside rectangles show only the rounded interval size.
    sr_label = format_genomic_length_label(sr_size_bp),
    
    lr_label = format_genomic_length_label(lrs_size_bp),
    
    is_single_exon_sr = !is.na(n_matched_exons) &
      n_matched_exons == 1L
  ) %>%
  mutate(
    # Axis extent is derived only from true genomic coordinates.
    true_rect_extent = pmax(
      abs(sr_xmin_true),
      abs(sr_xmax_true),
      abs(lr_xmin_true),
      abs(lr_xmax_true)
    ),
    
    # Exon-anchored uncertainty bounds are included in the panel extent.
    true_left_extent = if_else(
      left_bounded,
      abs(left_boundary_coordinate - cnv_midpoint),
      NA_real_
    ),
    
    true_right_extent = if_else(
      right_bounded,
      abs(right_boundary_coordinate - cnv_midpoint),
      NA_real_
    ),
    
    true_bounded_extent = pmax(
      true_rect_extent,
      coalesce(true_left_extent, 0),
      coalesce(true_right_extent, 0)
    ),
    
    any_unbounded = !left_bounded | !right_bounded,
    
    window_padding_multiplier = if_else(
      any_unbounded,
      1.35,
      1.05
    ),
    
    window_halfwidth = true_bounded_extent * window_padding_multiplier
  ) %>%
  mutate(
    chars_per_label = pmax(
      nchar(sr_label),
      nchar(lr_label)
    ),
    
    char_width_bp = (window_halfwidth * 2) / label_chars_per_window,
    
    min_display_width_bp = chars_per_label *
      char_width_bp *
      label_padding_multiplier,
    
    sr_true_width = sr_xmax_true - sr_xmin_true,
    lr_true_width = lr_xmax_true - lr_xmin_true,
    
    sr_true_centre = (sr_xmin_true + sr_xmax_true) / 2,
    
    # Widen only the short-read box for label legibility.
    sr_display_halfwidth = pmax(
      sr_true_width / 2,
      min_display_width_bp / 2
    ),
    
    sr_xmin = pmax(
      sr_true_centre - sr_display_halfwidth,
      -window_halfwidth
    ),
    
    sr_xmax = pmin(
      sr_true_centre + sr_display_halfwidth,
      window_halfwidth
    ),
    
    was_widened = is_single_exon_sr &
      (sr_display_halfwidth > (sr_true_width / 2)),
    
    # Long-read interval remains true to its genomic coordinates.
    lr_xmin = lr_xmin_true,
    lr_xmax = lr_xmax_true,
    
    sr_centre = (sr_xmin + sr_xmax) / 2,
    lr_centre = (lr_xmin + lr_xmax) / 2,
    
    # Critical correction: ticks use true adjacent-exon edge coordinates.
    left_tick_x = if_else(
      left_bounded,
      left_boundary_coordinate - cnv_midpoint,
      NA_real_
    ),
    
    right_tick_x = if_else(
      right_bounded,
      right_boundary_coordinate - cnv_midpoint,
      NA_real_
    ),
    
    # Uncertainty segments still connect to displayed SR-box edges.
    left_line_inner_x = sr_xmin,
    right_line_inner_x = sr_xmax
  ) %>%
  mutate(
    left_line_start = if_else(
      left_bounded,
      left_tick_x,
      -window_halfwidth
    ),
    
    right_line_end = if_else(
      right_bounded,
      right_tick_x,
      window_halfwidth
    ),
    
    left_uncertainty_label_x = (
      left_line_start + left_line_inner_x
    ) / 2,
    
    right_uncertainty_label_x = (
      right_line_inner_x + right_line_end
    ) / 2
  ) %>%
  select(
    -true_rect_extent,
    -true_left_extent,
    -true_right_extent,
    -true_bounded_extent
  ) %>%
  arrange(sample_id_trimmed) %>%
  mutate(
    panel_label = factor(
      panel_label,
      levels = unique(panel_label)
    )
  )

# --------------------- Coordinate audit -------------------------

# For bounded sides, the uncertainty tick coordinate must equal the
# corresponding adjacent exon boundary coordinate after centring.
uncertainty_audit <- cnv_data %>%
  transmute(
    sample_id_trimmed,
    Gene,
    left_bounded,
    left_boundary_coordinate,
    left_tick_x,
    left_tick_coordinate_check = left_tick_x + cnv_midpoint,
    left_tick_matches_exon_edge = !left_bounded |
      abs(left_tick_coordinate_check - left_boundary_coordinate) < 1e-8,
    right_bounded,
    right_boundary_coordinate,
    right_tick_x,
    right_tick_coordinate_check = right_tick_x + cnv_midpoint,
    right_tick_matches_exon_edge = !right_bounded |
      abs(right_tick_coordinate_check - right_boundary_coordinate) < 1e-8
  )

print(uncertainty_audit)

# Audit short-read widening.
print(
  cnv_data %>%
    select(
      sample_id_trimmed,
      Gene,
      n_matched_exons,
      is_single_exon_sr,
      sr_true_width,
      min_display_width_bp,
      was_widened,
      sr_display_halfwidth
    ) %>%
    mutate(
      widening_ratio = (2 * sr_display_halfwidth) / sr_true_width
    )
)

# ----------------- Transcript/exon/UTR plot data ----------------

gene_transcript_span <- exons %>%
  group_by(gene, chromosomename, strand) %>%
  summarise(
    transcript_start = min(exonchromstart),
    transcript_end = max(exonchromend),
    .groups = "drop"
  )

exon_track_data <- cnv_data %>%
  select(
    sample_id_trimmed,
    panel_label,
    Gene,
    lrs_chrom,
    cnv_midpoint,
    window_halfwidth
  ) %>%
  left_join(
    exons %>%
      select(
        gene,
        chromosomename,
        exon_rank,
        exonchromstart,
        exonchromend,
        strand
      ),
    by = c(
      "Gene" = "gene",
      "lrs_chrom" = "chromosomename"
    )
  ) %>%
  mutate(
    exon_xmin = pmax(exonchromstart - cnv_midpoint, -window_halfwidth),
    exon_xmax = pmin(exonchromend - cnv_midpoint, window_halfwidth)
  ) %>%
  filter(exon_xmax > exon_xmin)

utr_track_data <- exon_track_data %>%
  select(
    sample_id_trimmed,
    panel_label,
    Gene,
    lrs_chrom,
    cnv_midpoint,
    window_halfwidth,
    exon_rank,
    exon_xmin,
    exon_xmax
  ) %>%
  inner_join(
    utr_segments %>%
      select(
        gene,
        chromosomename,
        exon_rank,
        utr_start,
        utr_end,
        utr_type
      ),
    by = c(
      "Gene" = "gene",
      "lrs_chrom" = "chromosomename",
      "exon_rank" = "exon_rank"
    ),
    relationship = "many-to-many"
  ) %>%
  mutate(
    xmin = pmax(utr_start - cnv_midpoint, -window_halfwidth),
    xmax = pmin(utr_end - cnv_midpoint, window_halfwidth)
  ) %>%
  filter(xmax > xmin) %>%
  distinct(
    sample_id_trimmed,
    panel_label,
    Gene,
    lrs_chrom,
    exon_rank,
    utr_type,
    xmin,
    xmax,
    .keep_all = TRUE
  )

coding_exon_track_data <- exon_track_data %>%
  left_join(
    utr_track_data %>%
      group_by(
        sample_id_trimmed,
        panel_label,
        Gene,
        lrs_chrom,
        exon_rank
      ) %>%
      summarise(
        utr_left_edge = min(xmin),
        utr_right_edge = max(xmax),
        .groups = "drop"
      ),
    by = c(
      "sample_id_trimmed",
      "panel_label",
      "Gene",
      "lrs_chrom",
      "exon_rank"
    )
  ) %>%
  mutate(
    coding_xmin = if_else(
      !is.na(utr_left_edge) &
        abs(utr_left_edge - exon_xmin) < 1,
      utr_right_edge,
      exon_xmin
    ),
    
    coding_xmax = if_else(
      !is.na(utr_right_edge) &
        abs(utr_right_edge - exon_xmax) < 1,
      utr_left_edge,
      exon_xmax
    )
  ) %>%
  transmute(
    sample_id_trimmed,
    panel_label,
    Gene,
    lrs_chrom,
    exon_rank,
    xmin = coding_xmin,
    xmax = coding_xmax
  ) %>%
  filter(
    !is.na(xmin),
    !is.na(xmax),
    xmax > xmin
  )

# ----------------- Transcript backbone data ---------------------

transcript_data <- cnv_data %>%
  select(
    sample_id_trimmed,
    panel_label,
    Gene,
    lrs_chrom,
    cnv_midpoint,
    window_halfwidth
  ) %>%
  left_join(
    gene_transcript_span,
    by = c(
      "Gene" = "gene",
      "lrs_chrom" = "chromosomename"
    )
  ) %>%
  mutate(
    transcript_start_x = pmax(
      transcript_start - cnv_midpoint,
      -window_halfwidth
    ),
    
    transcript_end_x = pmin(
      transcript_end - cnv_midpoint,
      window_halfwidth
    )
  ) %>%
  filter(transcript_end_x > transcript_start_x)

# ----------------------- Track labels ---------------------------

transcript_y <- 0.30

# Track labels are positioned immediately to the right of each panel.
# Long-read and short-read label colours correspond to their interval tracks.
track_labels <- cnv_data %>%
  transmute(
    sample_id_trimmed,
    panel_label,
    label_x = window_halfwidth * 1.03,
    
    longread_label = "Long-read",
    longread_label_y = lr_mid_y,
    
    shortread_label = "Short-read",
    shortread_label_y = sr_mid_y,
    
    transcript_label = paste(
      Gene,
      unname(clinically_relevant_transcripts[Gene])
    ),
    transcript_label_y = transcript_y
  )

# --------------------- Plot coordinates -------------------------

# Long-read interval above short-read interval.
lr_y_min <- 1.55
lr_y_max <- 2.05

sr_y_min <- 0.85
sr_y_max <- 1.35

lr_mid_y <- (lr_y_min + lr_y_max) / 2
sr_mid_y <- (sr_y_min + sr_y_max) / 2

# Uncertainty ticks remain associated with the short-read call.
tick_y_min <- sr_y_min - 0.08
tick_y_max <- sr_y_max + 0.08

# Dotted guides terminate at the upper edge of exon rectangles.
sr_guide_y_min <- exon_track_y_max
sr_guide_y_max <- sr_y_min

sr_fill_colour <- "pink"
lr_fill_colour <- "lightblue"
exon_fill_colour <- "forestgreen"
utr_fill_colour <- "forestgreen"

# --------------------------- Plot -------------------------------

breakpoint_plot <- ggplot(cnv_data) +
  geom_blank(aes(x = -window_halfwidth)) +
  geom_blank(aes(x = window_halfwidth)) +
  
  # Left uncertainty range.
  geom_segment(
    data = cnv_data %>% filter(left_bounded),
    aes(
      x = left_line_start,
      xend = left_line_inner_x,
      y = sr_mid_y,
      yend = sr_mid_y
    ),
    linewidth = 1.2,
    colour = "grey45",
    lineend = "butt"
  ) +
  
  geom_segment(
    data = cnv_data %>% filter(!left_bounded),
    aes(
      x = left_line_inner_x,
      xend = left_line_start,
      y = sr_mid_y,
      yend = sr_mid_y
    ),
    linewidth = 1.2,
    colour = "grey45",
    arrow = arrow(
      type = "closed",
      length = grid::unit(3, "mm")
    )
  ) +
  
  # Right uncertainty range.
  geom_segment(
    data = cnv_data %>% filter(right_bounded),
    aes(
      x = right_line_inner_x,
      xend = right_line_end,
      y = sr_mid_y,
      yend = sr_mid_y
    ),
    linewidth = 1.2,
    colour = "grey45",
    lineend = "butt"
  ) +
  
  geom_segment(
    data = cnv_data %>% filter(!right_bounded),
    aes(
      x = right_line_inner_x,
      xend = right_line_end,
      y = sr_mid_y,
      yend = sr_mid_y
    ),
    linewidth = 1.2,
    colour = "grey45",
    arrow = arrow(
      type = "closed",
      length = grid::unit(3, "mm")
    )
  ) +
  
  # Bounded uncertainty ticks, now exactly exon-edge anchored.
  geom_segment(
    data = cnv_data %>% filter(left_bounded),
    aes(
      x = left_line_start,
      xend = left_line_start,
      y = tick_y_min,
      yend = tick_y_max
    ),
    linewidth = 1.0,
    colour = "grey30"
  ) +
  
  geom_segment(
    data = cnv_data %>% filter(right_bounded),
    aes(
      x = right_line_end,
      xend = right_line_end,
      y = tick_y_min,
      yend = tick_y_max
    ),
    linewidth = 1.0,
    colour = "grey30"
  ) +
  
  geom_text(
    data = cnv_data,
    aes(
      x = left_uncertainty_label_x,
      y = sr_mid_y + 0.10,
      label = left_uncertainty_label
    ),
    hjust = 0.5,
    vjust = 0,
    size = 2.8,
    colour = "grey20"
  ) +
  
  geom_text(
    data = cnv_data,
    aes(
      x = right_uncertainty_label_x,
      y = sr_mid_y + 0.10,
      label = right_uncertainty_label
    ),
    hjust = 0.5,
    vjust = 0,
    size = 2.8,
    colour = "grey20"
  ) +
  
  # Long-read interval.
  geom_rect(
    aes(
      xmin = lr_xmin,
      xmax = lr_xmax,
      ymin = lr_y_min,
      ymax = lr_y_max
    ),
    fill = lr_fill_colour,
    colour = "black",
    linewidth = 0.5
  ) +
  
  geom_text(
    aes(
      x = (lr_xmin + lr_xmax) / 2,
      y = lr_mid_y,
      label = lr_label
    ),
    size = 3.0,
    fontface = "bold",
    colour = "black"
  ) +
  
  # Short-read interval.
  geom_rect(
    aes(
      xmin = sr_xmin,
      xmax = sr_xmax,
      ymin = sr_y_min,
      ymax = sr_y_max
    ),
    fill = sr_fill_colour,
    colour = "black",
    linewidth = 0.5
  ) +
  
  # Dotted guides from displayed SR edges to true SR genomic edges.
  geom_segment(
    data = cnv_data,
    aes(
      x = sr_xmin,
      xend = sr_xmin_true,
      y = sr_guide_y_max,
      yend = sr_guide_y_min
    ),
    linewidth = 0.45,
    linetype = "dotted",
    colour = "grey35"
  ) +
  
  geom_segment(
    data = cnv_data,
    aes(
      x = sr_xmax,
      xend = sr_xmax_true,
      y = sr_guide_y_max,
      yend = sr_guide_y_min
    ),
    linewidth = 0.45,
    linetype = "dotted",
    colour = "grey35"
  ) +
  
  geom_text(
    aes(
      x = (sr_xmin + sr_xmax) / 2,
      y = sr_mid_y,
      label = sr_label
    ),
    size = 3.0,
    fontface = "bold",
    colour = "black"
  ) +
  
  # Transcript backbone.
  geom_segment(
    data = transcript_data,
    aes(
      x = transcript_start_x,
      xend = transcript_end_x,
      y = transcript_y,
      yend = transcript_y
    ),
    linewidth = 0.4,
    colour = "darkgreen"
  ) +
  
  # Coding exon regions.
  geom_rect(
    data = coding_exon_track_data,
    aes(
      xmin = xmin,
      xmax = xmax,
      ymin = exon_track_y_min,
      ymax = exon_track_y_max
    ),
    fill = exon_fill_colour,
    colour = "black",
    linewidth = 0.2
  ) +
  
  # UTR regions.
  geom_rect(
    data = utr_track_data,
    aes(
      xmin = xmin,
      xmax = xmax,
      ymin = utr_track_y_min,
      ymax = utr_track_y_max
    ),
    fill = utr_fill_colour,
    colour = "black",
    linewidth = 0.2
  ) +
  
  # Per-panel track labels.
  geom_text(
    data = track_labels,
    aes(
      x = label_x,
      y = longread_label_y,
      label = longread_label
    ),
    hjust = 0,
    vjust = 0.5,
    fontface = "bold",
    size = 3.0,
    colour = "lightblue"
  ) +
  
  geom_text(
    data = track_labels,
    aes(
      x = label_x,
      y = shortread_label_y,
      label = shortread_label
    ),
    hjust = 0,
    vjust = 0.5,
    fontface = "bold",
    size = 3.0,
    colour = "pink"
  ) +
  
  geom_text(
    data = track_labels,
    aes(
      x = label_x,
      y = transcript_label_y,
      label = transcript_label
    ),
    hjust = 0,
    vjust = 0.5,
    size = 3.0,
    colour = "forestgreen"
  ) +
  
  facet_wrap(
    ~panel_label,
    ncol = 1,
    scales = "free_x",
    strip.position = "left"
  ) +
  
  scale_y_continuous(
    limits = c(0.15, 2.15)
  ) +
  
  scale_x_continuous(
    expand = expansion(mult = c(0.08, 0.08))
  ) +
  
  coord_cartesian(clip = "off") +
  
  labs(
    title = "Uncertainty of short-read CNV interval vs long-read CNV interval",
    subtitle = paste(
      "Terminal-exon UTRs are shown as thinner boxes. Dotted guides link",
      "short-read interval edges to their caller-assigned genomic breakpoint coordinates."
    ),
    x = NULL,
    y = NULL
  ) +
  
  theme_classic(base_size = 11) +
  
  theme(
    legend.position = "none",
    
    axis.line.x = element_blank(),
    axis.line.y = element_blank(),
    
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_blank(),
    
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    
    panel.border = element_rect(
      colour = "grey55",
      fill = NA,
      linewidth = 0.4
    ),
    
    strip.text.y.left = element_text(
      angle = 0,
      hjust = 1,
      size = 8.5,
      lineheight = 0.9
    ),
    
    strip.background = element_rect(
      fill = "grey97",
      colour = "grey55",
      linewidth = 0.4
    ),
    
    panel.spacing = grid::unit(1.2, "lines"),
    
    plot.title = element_text(
      face = "bold",
      size = 14
    ),
    
    plot.subtitle = element_text(
      size = 9,
      colour = "grey40"
    ),
    
    plot.margin = margin(10, 50, 8, 15)
  )

breakpoint_plot

# ---------------- Page geometry for landscape A4 ----------------

a4_landscape_width_in <- 11.72
a4_landscape_height_in <- 8.27

page_margin_lr_in <- 0.5
page_margin_tb_in <- 0.4

caption_font_pt <- 12
caption_lines <- 3
caption_line_height <- 1.15
caption_gap_above_in <- 0.15

caption_height_in <- (
  caption_font_pt *
    caption_line_height *
    caption_lines
) / 72 + caption_gap_above_in

usable_width_in <- a4_landscape_width_in - page_margin_lr_in

usable_height_in <- a4_landscape_height_in -
  page_margin_tb_in -
  caption_height_in

desired_height_in <- 1.9 *
  length(unique(cnv_data$panel_label)) +
  1.5

final_height_in <- min(
  desired_height_in,
  usable_height_in
)

final_width_in <- usable_width_in

if (desired_height_in > usable_height_in) {
  warning(
    "Figure height (",
    round(desired_height_in, 2),
    " in) exceeds usable landscape A4 space (",
    round(usable_height_in, 2),
    " in)."
  )
}

# ggsave(
#   filename = paste0(output_prefix, "_plot.png"),
#   plot = breakpoint_plot,
#   width = final_width_in,
#   height = final_height_in,
#   units = "in",
#   dpi = 600
# )
