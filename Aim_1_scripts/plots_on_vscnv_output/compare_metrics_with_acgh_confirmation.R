## This plot is for Section 3 of my thesis where it outputs the 
## naive CNVs that are within the simulated range
## AND show which CNVs are selected for validation!

library(ggplot2)
library(plotly)
library(dplyr)
library(ggpubr)
library(pROC)

real_data <- sample_cnv_calls_combined

validated_res <- read.csv("/Users/zlaw0001/Library/CloudStorage/OneDrive-MonashUniversity/Joshua_PhD_Project/Aim_1/Validations/acgh_validation_status_files/acgh_round_one_res_for_candidate_cnvs_dan.csv")

## omitting the maybe, likely and "a" because they are CNV calls that Dan cannot be certain if
## concordant or discordant.
# validated_res <- validated_res %>% filter(!Verified_20250730 %in% c("MAYBE", "Likely", "a"))

validated_res$Man..cur..Chr <- gsub("^chr", "", validated_res$Man..cur..Chr)  # remove 'chr' prefix
validated_res$Man..cur..Chr[validated_res$Man..cur..Chr == ""] <- NA            # convert empty strings to NA

validated_res$Var.Seq.Chr <- gsub("^chr", "", validated_res$Var.Seq.Chr)  # remove 'chr' prefix
validated_res$Var.Seq.Chr[validated_res$Var.Seq.Chr == ""] <- NA      

validated_res[] <- lapply(validated_res, function(x) {
  if (is.character(x)) x[x == ""] <- NA
  return(x)
})

validated_res <- validated_res %>%
  mutate(VS_Call_chr = as.character(VS_Call_chr))

real_data_with_acgh <- real_data %>%
  left_join(
    validated_res %>%
      select(sample,
             VS_Call_chr,
             VS_Call_start,
             VS_Call_end,
             Array.Global.Display.name,
             Verified_20250730),
    by = c("sample" = "sample",
           "VS_Call_chr" = "VS_Call_chr",
           "VS_Call_start" = "VS_Call_start",
           "VS_Call_end" = "VS_Call_end")
  )


real_data_with_acgh <- real_data_with_acgh %>%
  mutate(passed_score_filters = ifelse((cn == 3 & 
                                          (Avg.Read.Ratio.of.CNV > min_clinvar_dup_rr_diff) & 
                                          (Avg.Z.Score.of.CNV > min_clinvar_dup_z_score)) |
                                         (cn == 1  & 
                                            (Avg.Read.Ratio.of.CNV < max_clinvar_del_rr_diff) & 
                                            (Avg.Z.Score.of.CNV < max_clinvar_del_z_score)), TRUE, FALSE)) %>%
  mutate(fail_percent_diff = ifelse(Percent.Difference > 19, TRUE, FALSE))


# filter out CNVs from outlier samples

real_data_with_acgh <- real_data_with_acgh %>%
  filter(outlier_sample == "NO") %>%
  filter(fail_percent_diff == FALSE)
  # filter(!Verified_20250730 %in% c("MAYBE", "Likely", "a", NA))



#### READ RATIO vs Z-score ####

#### SELECTED and ran ACGH ####
real_data_with_acgh_selected_vs_unselected <- real_data_with_acgh %>%
  mutate(color_group = case_when(
    # sent_for_acgh == TRUE ~ "Selected for aCGH Validation",
    !is.na(Verified_20250730) ~ "Sent for aCGH validation",
    is.na(Verified_20250730) ~ "Not sent for aCGH validation",
    TRUE ~ "gray"
  ))


real_data_with_acgh_rr_vs_z_score <- ggplot(real_data_with_acgh_selected_vs_unselected,
                                            aes(x = Avg.Read.Ratio.of.CNV,
                                                y = Avg.Z.Score.of.CNV,
                                                color = color_group, text = text)) +
  geom_point(size = 3, alpha = 0.6) +
  
  # Vertical lines for max/min clinvar Read Ratios
  geom_vline(xintercept = max_clinvar_del_rr_diff, linetype = "dashed", color = "black") +
  geom_vline(xintercept = min_clinvar_dup_rr_diff, linetype = "dashed", color = "black") +
  
  # Horizontal lines for max/min clinvar Z-scores
  geom_hline(yintercept = max_clinvar_del_z_score, linetype = "dashed", color = "black") +
  geom_hline(yintercept = min_clinvar_dup_z_score, linetype = "dashed", color = "black") +
  
  # Annotations for Read Ratio thresholds (top of plot)
  annotate("text", x = max_clinvar_del_rr_diff - 0.3,
           y = max(real_data_without_outliers$Avg.Z.Score.of.CNV, na.rm = TRUE),
           label = paste0("Simulated deletion read ratio maximum = ", round(max_clinvar_del_rr_diff, 2)),
           color = "black", vjust = -0.75, size = 4.5) +
  annotate("text", x = min_clinvar_dup_rr_diff + 0.3,
           y = max(real_data_without_outliers$Avg.Z.Score.of.CNV, na.rm = TRUE),
           label = paste0("Simulated duplication read ratio minimum = ", round(min_clinvar_dup_rr_diff, 2)),
           color = "black", vjust = -0.75, size = 4.5) +
  
  # Annotations for Z-score thresholds (sides of plot)
  annotate("text", x = min(real_data_without_outliers$Avg.Read.Ratio.of.CNV, na.rm = TRUE),
           y = max_clinvar_del_z_score,
           label = paste0("Simulated deletion z-score maximum = ", round(max_clinvar_del_z_score, 2)),
           color = "black", vjust = -1, hjust = 0, size = 4.5) +
  annotate("text", x = min(real_data_without_outliers$Avg.Read.Ratio.of.CNV, na.rm = TRUE),
           y = min_clinvar_dup_z_score,
           label = paste0("Simulated duplication z-score minimum = ", round(min_clinvar_dup_z_score, 2)),
           color = "black", vjust = 1.5, hjust = 0, size = 4.5) +
  
  # Labels and styling
  labs(
    x = "Read Ratio",
    y = "Z-Score",
    color = "CNV Type: "
  ) +
  scale_color_manual(values = c(
    'Sent for aCGH validation' = 'purple',
    'Not sent for aCGH validation' = 'pink'
  )) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.title = element_text(size = 13, face = "bold"),
    legend.text = element_text(size = 12),
    legend.key.size = unit(0.8, "lines"),
    legend.spacing.x = unit(0.4, 'cm'),
    
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    axis.title.x = element_text(size = 14, face = "bold"),
    axis.title.y = element_text(size = 14, face = "bold"),
    axis.text = element_text(size = 12)
  )

real_data_with_acgh_rr_vs_z_score

# Convert to ggplotly
real_data_with_acgh_rr_vs_z_score_plotly <- ggplotly(real_data_with_acgh_rr_vs_z_score, tooltip = "text")
real_data_with_acgh_rr_vs_z_score_plotly



#### WITH ACGH CONFIRMATION ####
real_data_with_acgh_with_acgh_confirmation <- real_data_with_acgh %>%
  filter(!is.na(Verified_20250730))

real_data_with_acgh_with_acgh_confirmation <- real_data_with_acgh_with_acgh_confirmation %>%
  mutate(color_group = case_when(
    # sent_for_acgh == TRUE ~ "Selected for aCGH Validation",
    Verified_20250730 == "YES" ~ "Supported by aCGH",
    !Verified_20250730 %in%  c("MAYBE", "Likely", "a", "YES", NA) ~ "NOT supported by aCGH",
    Verified_20250730 %in% c("MAYBE", "Likely", "a") ~ "Inconclusive",
    TRUE ~ "gray"
  ))


real_data_with_acgh_rr_vs_z_score <- ggplot(real_data_with_acgh_with_acgh_confirmation,
                                       aes(x = Avg.Read.Ratio.of.CNV,
                                           y = Avg.Z.Score.of.CNV,
                                           color = color_group, text = text)) +
  geom_point(size = 3, alpha = 0.6) +
  
  # Vertical lines for max/min clinvar Read Ratios
  geom_vline(xintercept = max_clinvar_del_rr_diff, linetype = "dashed", color = "black") +
  geom_vline(xintercept = min_clinvar_dup_rr_diff, linetype = "dashed", color = "black") +
  
  # Horizontal lines for max/min clinvar Z-scores
  geom_hline(yintercept = max_clinvar_del_z_score, linetype = "dashed", color = "black") +
  geom_hline(yintercept = min_clinvar_dup_z_score, linetype = "dashed", color = "black") +
  
  # Annotations for Read Ratio thresholds (top of plot)
  annotate("text", x = max_clinvar_del_rr_diff - 0.3,
           y = max(real_data_without_outliers$Avg.Z.Score.of.CNV, na.rm = TRUE),
           label = paste0("Simulated deletion read ratio maximum = ", round(max_clinvar_del_rr_diff, 2)),
           color = "black", vjust = -0.75, size = 4.5) +
  annotate("text", x = min_clinvar_dup_rr_diff + 0.3,
           y = max(real_data_without_outliers$Avg.Z.Score.of.CNV, na.rm = TRUE),
           label = paste0("Simulated duplication read ratio minimum = ", round(min_clinvar_dup_rr_diff, 2)),
           color = "black", vjust = -0.75, size = 4.5) +
  
  # Annotations for Z-score thresholds (sides of plot)
  annotate("text", x = min(real_data_without_outliers$Avg.Read.Ratio.of.CNV, na.rm = TRUE),
           y = max_clinvar_del_z_score,
           label = paste0("Simulated deletion z-score maximum = ", round(max_clinvar_del_z_score, 2)),
           color = "black", vjust = -1, hjust = 0, size = 4.5) +
  annotate("text", x = min(real_data_without_outliers$Avg.Read.Ratio.of.CNV, na.rm = TRUE),
           y = min_clinvar_dup_z_score,
           label = paste0("Simulated duplication z-score minimum = ", round(min_clinvar_dup_z_score, 2)),
           color = "black", vjust = 1.5, hjust = 0, size = 4.5) +
  
  # Labels and styling
  labs(
    x = "Read Ratio",
    y = "Z-Score",
    color = "CNV Type: "
  ) +
  scale_color_manual(values = c(
    'Supported by aCGH' = '#2196F3',
    'NOT supported by aCGH' = '#F44336',
    'Inconclusive' = 'grey'
  )) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.title = element_text(size = 13, face = "bold"),
    legend.text = element_text(size = 12),
    legend.key.size = unit(0.8, "lines"),
    legend.spacing.x = unit(0.4, 'cm'),
    
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    axis.title.x = element_text(size = 14, face = "bold"),
    axis.title.y = element_text(size = 14, face = "bold"),
    axis.text = element_text(size = 12)
  )

real_data_with_acgh_rr_vs_z_score

# Convert to ggplotly
real_data_with_acgh_rr_vs_z_score_plotly <- ggplotly(real_data_with_acgh_rr_vs_z_score, tooltip = "text")
real_data_with_acgh_rr_vs_z_score_plotly



#### Selection of CNVs for aCGH validation and their results #####

#### DELETIONS - Z-score ranks with top-100 tracking column ####

deletion_z_ranks <- real_data_with_acgh %>%
  filter(Avg.Read.Ratio.of.CNV < 1.0) %>%
  mutate(
    z_score_rank   = rank(Avg.Z.Score.of.CNV, ties.method = "min"),
    top_100_status = ifelse(z_score_rank <= 100, "Top-100", "Outside top-100"),
    acgh_outcome   = case_when(
      Verified_20250730 == "YES"                              ~ "Confirmed",
      Verified_20250730 %in% c("MAYBE", "Likely", "a")         ~ "Inconclusive",
      !is.na(Verified_20250730)                                ~ "Unsupported",
      TRUE                                                      ~ "Not tested"
    )
  )

deletion_summary <- deletion_z_ranks %>%
  filter(acgh_outcome != "Not tested") %>%
  group_by(top_100_status, acgh_outcome) %>%
  summarise(n = n(), .groups = "drop") %>%
  tidyr::pivot_wider(names_from = acgh_outcome, values_from = n, values_fill = 0) %>%
  mutate(Total_tested = Confirmed + Unsupported + Inconclusive)

deletion_sim_fail_summary <- deletion_z_ranks %>%
  filter(acgh_outcome != "Not tested", passed_score_filters != TRUE) %>%
  group_by(top_100_status, acgh_outcome) %>%
  summarise(n = n(), .groups = "drop") %>%
  tidyr::pivot_wider(names_from = acgh_outcome, values_from = n, values_fill = 0)

print(deletion_summary)
print(deletion_sim_fail_summary)


#### DUPLICATIONS - Z-score ranks with top-100 tracking column ####

duplication_z_ranks <- real_data_with_acgh %>%
  filter(Avg.Read.Ratio.of.CNV > 1.0) %>%
  mutate(
    z_score_rank   = rank(-Avg.Z.Score.of.CNV, ties.method = "min"),
    top_100_status = ifelse(z_score_rank <= 100, "Top-100", "Outside top-100"),
    acgh_outcome   = case_when(
      Verified_20250730 == "YES"                              ~ "Confirmed",
      Verified_20250730 %in% c("MAYBE", "Likely", "a")         ~ "Inconclusive",
      !is.na(Verified_20250730)                                ~ "Unsupported",
      TRUE                                                      ~ "Not tested"
    )
  )

duplication_summary <- duplication_z_ranks %>%
  filter(acgh_outcome != "Not tested") %>%
  group_by(top_100_status, acgh_outcome) %>%
  summarise(n = n(), .groups = "drop") %>%
  tidyr::pivot_wider(names_from = acgh_outcome, values_from = n, values_fill = 0) %>%
  mutate(Total_tested = Confirmed + Unsupported + Inconclusive)

duplication_sim_fail_summary <- duplication_z_ranks %>%
  filter(acgh_outcome != "Not tested", passed_score_filters != TRUE) %>%
  group_by(top_100_status, acgh_outcome) %>%
  summarise(n = n(), .groups = "drop") %>%
  tidyr::pivot_wider(names_from = acgh_outcome, values_from = n, values_fill = 0)

print(duplication_summary)
print(duplication_sim_fail_summary)


#### Combined table across both CNV types, for reporting in Section 3 ####

combined_z_ranks <- bind_rows(
  deletion_z_ranks    %>% mutate(cnv_type = "Deletion"),
  duplication_z_ranks %>% mutate(cnv_type = "Duplication")
)

combined_summary <- combined_z_ranks %>%
  filter(acgh_outcome != "Not tested") %>%
  group_by(cnv_type, top_100_status, acgh_outcome) %>%
  summarise(n = n(), .groups = "drop") %>%
  tidyr::pivot_wider(names_from = acgh_outcome, values_from = n, values_fill = 0) %>%
  mutate(Total_tested = Confirmed + Unsupported + Inconclusive)

print(combined_summary)


#### DELETIONS - Z-score and Read Ratio ranges (aCGH validated only) ####

deletion_acgh_ranges <- real_data_with_acgh %>%
  filter(Avg.Read.Ratio.of.CNV < 1.0, !is.na(Verified_20250730)) %>%
  summarise(
    min_z_score = min(Avg.Z.Score.of.CNV, na.rm = TRUE),
    max_z_score = max(Avg.Z.Score.of.CNV, na.rm = TRUE),
    min_read_ratio = min(Avg.Read.Ratio.of.CNV, na.rm = TRUE),
    max_read_ratio = max(Avg.Read.Ratio.of.CNV, na.rm = TRUE)
  )

print(paste0("Deletions (aCGH validated) Z-score range: ",
             round(deletion_acgh_ranges$min_z_score, 4), " to ",
             round(deletion_acgh_ranges$max_z_score, 4)))

print(paste0("Deletions (aCGH validated) Read Ratio range: ",
             round(deletion_acgh_ranges$min_read_ratio, 4), " to ",
             round(deletion_acgh_ranges$max_read_ratio, 4)))



#### DUPLICATIONS - Z-score and Read Ratio ranges (aCGH validated only) ####

duplication_acgh_ranges <- real_data_with_acgh %>%
  filter(Avg.Read.Ratio.of.CNV > 1.0, !is.na(Verified_20250730)) %>%
  summarise(
    min_z_score = min(Avg.Z.Score.of.CNV, na.rm = TRUE),
    max_z_score = max(Avg.Z.Score.of.CNV, na.rm = TRUE),
    min_read_ratio = min(Avg.Read.Ratio.of.CNV, na.rm = TRUE),
    max_read_ratio = max(Avg.Read.Ratio.of.CNV, na.rm = TRUE)
  )

print(paste0("Duplications (aCGH validated) Z-score range: ",
             round(duplication_acgh_ranges$min_z_score, 4), " to ",
             round(duplication_acgh_ranges$max_z_score, 4)))

print(paste0("Duplications (aCGH validated) Read Ratio range: ",
             round(duplication_acgh_ranges$min_read_ratio, 4), " to ",
             round(duplication_acgh_ranges$max_read_ratio, 4)))




#### Median summary: Z-score and Read Ratio by CNV class and aCGH outcome ####

median_summary <- auroc_data %>%
  group_by(cnv_class, color_group) %>%
  summarise(
    n               = n(),
    median_z_score  = median(Avg.Z.Score.of.CNV, na.rm = TRUE),
    median_read_ratio = median(Avg.Read.Ratio.of.CNV, na.rm = TRUE),
    .groups = "drop"
  )

print(median_summary)


######## AUROC with z-score and read ratio for aCGH confirmed and unsupported calls ########

## Restrict to calls with a definitive aCGH outcome (exclude inconclusive/NA)
auroc_data <- real_data_with_acgh %>%
  filter(color_group %in% c("Supported by aCGH", "NOT supported by aCGH")) %>%
  mutate(cnv_class = ifelse(Avg.Read.Ratio.of.CNV < 1, "Deletion", "Duplication"),
         outcome = ifelse(color_group == "Supported by aCGH", 1, 0))

#### DELETIONS ####
del_data <- auroc_data %>% filter(cnv_class == "Deletion")

roc_del_z  <- roc(del_data$outcome, del_data$Avg.Z.Score.of.CNV,      direction = ">")
roc_del_rr <- roc(del_data$outcome, del_data$Avg.Read.Ratio.of.CNV,   direction = ">")

auc_del_z  <- auc(roc_del_z)
auc_del_rr <- auc(roc_del_rr)

#### DUPLICATIONS ####
dup_data <- auroc_data %>% filter(cnv_class == "Duplication")

roc_dup_z  <- roc(dup_data$outcome, dup_data$Avg.Z.Score.of.CNV,      direction = "<")
roc_dup_rr <- roc(dup_data$outcome, dup_data$Avg.Read.Ratio.of.CNV,   direction = "<")

auc_dup_z  <- auc(roc_dup_z)
auc_dup_rr <- auc(roc_dup_rr)

#### Summary table ####
auroc_summary <- data.frame(
  CNV_Type = c("Deletion", "Deletion", "Duplication", "Duplication"),
  Metric   = c("Z-score", "Read ratio", "Z-score", "Read ratio"),
  AUROC    = round(c(auc_del_z, auc_del_rr, auc_dup_z, auc_dup_rr), 3)
)

print(auroc_summary)

#### Optional: plot ROC curves ####
plot(roc_del_z,  main = "Deletion Z-score ROC", col = "blue")
plot(roc_del_rr, main = "Deletion Read Ratio ROC", col = "red")
plot(roc_dup_z,  main = "Duplication Z-score ROC", col = "blue")
plot(roc_dup_rr, main = "Duplication Read Ratio ROC", col = "red")





