library(dplyr)
library(GenomicRanges)
library(tidyr)

##### THIS SCRIPT IS OBSOLETE NOW THAT I IMPLEMENTED assign_cyto_overlap FUNCTION ##########
##### ALSO THE CONCLUSION WE HAVE FOR DECON CLEARCNV AND GATK IS THAT NONE OF THESE TOOLS HAD ANY NON-VSCNV CALLS THAT OVERLAPPED
#### WITH ANY ACGH CALL FROM THE CYTO TABLES. SO WE ARE NOT MISSING ANY CNV THAT IS POTENTIALLY CALLED BY ACGH BUT NOT A CANDIDATE BY
### VS-CNV WHICH DAN IS REFERRING TO.

##### the goal of this script is to give Dan a list of CNVs that are NOT CANDIDATE VS-CNV CALLS AND OVERLAP WITH HIS CYTO
##### MASTER TABLE!
## Okay i've checked based on sample-gene pairs, and there is only one overlap, but that
## pair has different coordinates to what the aCGH sample-gene pair has, so there is no
## non-varseq calls from "decon", "clearCNV", "gatk" that overlap with any acgh calls

gene_whitelist <- c("BRCA1", "BRCA2", "PALB2", "MLH1", "MSH2", "MSH6", "APOB", "PCSK9", "LDLR")

# Prepare cyto table
cyto_master_table <- read.csv(
  "/path/to/project/Aim_1/Validations/stitching_dan_cyto_files/cyto_master_table.tsv", 
  sep = "\t"
) %>%
  mutate(
    sample_id_trimmed = gsub("^.*?_(.*?)_.*$", "\\1", sample_id),
    Chr = sub("^chr", "", Chr)
  ) %>%
  separate_rows(Gene.Names, sep = ",") %>%          # split genes into separate rows
  mutate(Gene.Names = trimws(Gene.Names)) %>%       # trim whitespace
  filter(Gene.Names %in% gene_whitelist)            # keep only whitelist genes


# Step 1: Convert cytogenetic CNVs to GRanges
cyto_gr <- GRanges(
  seqnames = cyto_master_table$Chr,
  ranges = IRanges(start = cyto_master_table$Start, end = cyto_master_table$Stop)
)

# Step 2: Combine all tool calls
all_calls <- bind_rows(lapply(names(harmonised_list_with_nonvscnv), function(tool_name) {
  if (tool_name %in% c("decon", "clearCNV", "gatk")) {
    df <- harmonised_list_with_nonvscnv[[tool_name]]$data
    df %>%
      filter(non_vscnv == "yes") %>%
      mutate(
        Tool = tool_name,
        sample_id_trimmed = sub(".*-(DNS\\d+)_S\\d+.*", "\\1", Sample)
      ) %>%
      select(Tool, Sample, sample_id_trimmed, Gene, Chromosome, Start, End, CNV_Type, everything())
  }
}))

decon_with_acgh$Tool <- "decon"
clearcnv_with_acgh$Tool <- "clearcnv"
gatk_with_acgh$Tool <- "gatk"

all_calls <- dplyr::bind_rows(
  decon_with_acgh,
  clearcnv_with_acgh,
  gatk_with_acgh
)

# Step 3: Filter tool calls to only sample-gene pairs present in cyto_master_table
all_calls_filtered <- all_calls %>%
  semi_join(
    cyto_master_table %>% 
      select(sample_id_trimmed, Gene = Gene.Names),
    by = c("sample_id_trimmed", "Gene")
  )

# Step 4: Convert tool calls to GRanges
calls_gr <- GRanges(
  seqnames = all_calls_filtered$Chromosome,
  ranges = IRanges(start = all_calls_filtered$Start, end = all_calls_filtered$End)
)

# Step 5: Subset cytogenetic GRanges to matching sample-gene pairs
cyto_gr_filtered <- cyto_gr[with(cyto_master_table, paste(sample_id_trimmed, Gene.Names)) %in% 
                              paste(all_calls_filtered$sample_id_trimmed, all_calls_filtered$Gene)]

# Step 6: Find overlaps
hits <- findOverlaps(calls_gr, cyto_gr_filtered)

# Step 7: Extract overlapping calls and annotate
overlapping_calls <- all_calls_filtered[queryHits(hits), ] %>%
  mutate(cyto_overlap = TRUE) %>%
  distinct(Tool, Sample, Chromosome, Start, End, CNV_Type, .keep_all = TRUE)

# Result
overlapping_calls