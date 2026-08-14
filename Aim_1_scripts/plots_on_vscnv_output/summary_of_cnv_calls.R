# Load necessary libraries
library(dplyr)
library(ggplot2)

# Number of Samples called per CNV
cnv_counts <- sample_cnv_calls_combined %>%
  group_by(coordinates, Estimated.CN.of.CNV) %>%   # Group by unique CNV calls
  summarise(
    Gene.Names = first(Gene.Names),  # Keep the first occurrence of Gene.Names
    PLP = first(PLP),                # Keep the first occurrence of PLP
    exon_type = first(exon_type),    # Keep the first occurrence of exon_type
    Number_of_Samples = n()          # Count the number of samples per CNV
  ) %>%
  arrange(desc(Number_of_Samples))   # Arrange by 'Number_of_Samples' in descending order

# View(cnv_counts)
# write.csv(cnv_counts, "summary_analyses_cnv_calls/cnv_counts.csv")

# histogram for cnv_counts
ggplot(cnv_counts, aes(x = Number_of_Samples)) +
  geom_histogram(binwidth = 1, fill = "blue", color = "black") + # Adjust binwidth as needed
  labs(title = "Histogram of Number of Samples CNV is called in",
       x = "Number of Samples",
       y = "Frequency") +
  theme_minimal()

# Number of CNVs called per gene

# Filter the data for CN=1 and CN=3 CNVs
cnv_per_gene <- sample_cnv_calls_combined %>%
  filter(Estimated.CN.of.CNV %in% c(1, 3)) %>%
  group_by(Gene.Names, Estimated.CN.of.CNV) %>%
  summarize(count = n()) %>%
  arrange(Gene.Names, Estimated.CN.of.CNV)

# View(cnv_per_gene)
# write.csv(cnv_per_gene, "summary_analyses_cnv_calls/cnv_per_gene.csv")
