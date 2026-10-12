library(dplyr)
library(ggplot2)

## Make sure vscnv_reorder_perrun is created from overlap_calls_from_tools_and_vscnv.R
# validated_res <- read.csv("/path/to/project/Aim_1/CNV_calling/cnvbenchmarker/CNVbenchmarkeR2_test/examples/final_validated_run_val_samples.tsv", sep = "\t")
validated_res <- read.csv("/path/to/project/Aim_1/Validations/acgh_validation_status_files/acgh_round_one_confirmed_cnvs_dan.csv")

validated_res$Man..cur..Chr <- gsub("^chr", "", validated_res$Man..cur..Chr)  # remove 'chr' prefix
validated_res$Man..cur..Chr[validated_res$Man..cur..Chr == ""] <- NA            # convert empty strings to NA

validated_res[] <- lapply(validated_res, function(x) {
  if (is.character(x)) x[x == ""] <- NA
  return(x)
})

validated_res_clean <- validated_res %>%
  rename(
    Sample = sample,
    Chromosome = VS_Call_chr,
    Start = VS_Call_start,
    End = VS_Call_end,
    Gene = Gene.Names,
    CNV_Type = Type.of.CNV
  ) %>%
  arrange(adjusted_p_value) %>% 
  mutate(Sample = sub("\\.hq\\.sorted\\.marked.*", "", Sample),
         Chromosome = gsub("^chr", "", Chromosome),
         CNV_Type = recode(CNV_Type,
                           "Duplicate" = "duplication",
                           "Het Deletion" = "deletion")) %>%
  filter(!is.na(Gene)) %>% ## Filter out any CNVs in ATM ROI
  dplyr::select(Sample, Chromosome, Start, End, Gene, CNV_Type, everything())

concordant <- validated_res_clean %>% filter(Verified_20250730 == "YES")
discordant <- validated_res_clean %>% filter(Verified_20250730 != "YES")

# concordant <- concordant %>%
#   mutate(depth_over_size_ratio = Avg.Target.Mean.Depth.of.CNV/Size.of.CNV)
# 
# discordant <- discordant %>%
#   mutate(depth_over_size_ratio = Avg.Target.Mean.Depth.of.CNV/Size.of.CNV)

single_exon_concordant <- concordant %>%
  filter(Number.of.Exons.by.VarSeq.call == 1)

single_exon_discordant <- discordant %>%
  filter(Number.of.Exons.by.VarSeq.call == 1)


# Add a label to each dataset
concordant <- concordant %>%
  mutate(Status = "Concordant")

discordant <- discordant %>%
  mutate(Status = "Discordant")

# Combine both
combined <- bind_rows(concordant, discordant)


## Generate a jitter plot of “number of samples with identical CNV call” for both concordant and discordant calls

ggplot(combined, aes(x = Status, y = Number_of_Samples, color = Status)) +
  geom_jitter(width = 0.2, alpha = 0.6, size = 2) +
  theme_minimal(base_size = 14) +
  labs(
    title = "Comparison of Number of Samples with identical CNV call in Concordant vs Discordant CNVs",
    x = "",
    y = "Number of Samples"
  ) +
  theme(legend.position = "none")



## Generate a jitter plot of “Avg.Target.Mean.Depth.of.CNV” for both concordant and discordant calls.

# Plot jitter comparison
ggplot(combined, aes(x = Status, y = Avg.Target.Mean.Depth.of.CNV, color = Status)) +
  geom_jitter(width = 0.2, alpha = 0.6, size = 2) +
  theme_minimal(base_size = 14) +
  labs(
    title = "Comparison of Average Target Mean Depth of CNV between Concordant and Discordant CNVs",
    x = "",
    y = "Average Target Mean Depth of CNV (Calculated by VS-CNV)"
  ) +
  theme(legend.position = "none")







