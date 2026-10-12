library(dplyr)
library(ggplot2)
library(tidyr)
library(IRanges)  # for overlap detection
library(stringr)
library(GenomicRanges)
library(patchwork)

# 1️⃣ Match factor order for both plots
zscore_combined$Tool <- factor(zscore_combined$Tool, levels = tool_order)


# Hide x labels for concordance
concordance_plot_adj <- concordance_plot +
  theme(
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    axis.ticks.x = element_blank()
  )

# Show x labels for jitter
combined_plot_adj <- violin_plot +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
    axis.title.x = element_blank(),
    panel.grid.major.x = element_blank()
  )

# Stack vertically
final_plot <- concordance_plot_adj / combined_plot_adj +
  plot_layout(heights = c(2, 1))

print(final_plot)

ggsave(
  filename = file.path(dir, "summary/combined_violet_concordance_hgsa.png"),
  plot = final_plot,
  width = 6.4*3,
  height = 3.6*3,
  dpi = 300,
  bg = "white"
)
