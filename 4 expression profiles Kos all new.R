#analyzing the expression of different (FunCs) in microbial genomes over time
# Load necessary libraries
library(tidyverse)

# Read in sample list, representative bins, and annotation
samples <- read.table("sample_list.txt")
regelist <- read.table("representative_bins.filtered.tax.tsv")
anno <- readRDS("Reges_annotations_genesfuncstax.RDS") %>%
  select(binID, FunC_ID, order, family, genus, species) %>%
  distinct() %>%
  mutate(FunC_ID = as.factor(FunC_ID))

# Read in expression data
expr <- readRDS("GeneLevel_Expressiondata_combined_TimeSeries.RDS")

# Subset data to only include certain bins and certain dates
expr <- expr %>%
  filter(binID %in% regelist$V1) %>%
  filter(date > "2011-01-25") %>%
  select(GeneID, binID, KO, date, MT_depth_scaledmappedreads, MG_depth_scaledmappedreads, MP_spectr_counts)

# Remove missing values and create expression profiles based on KO, ReGe
profiles <- expr %>%
  filter(!is.na(KO)) %>%
  group_by(binID, KO, date) %>%
  summarise(MT_mean = mean(MT_depth_scaledmappedreads),
            MG_mean = mean(MG_depth_scaledmappedreads),
            MP_sum = sum(MP_spectr_counts, na.rm=T),
            MT_max = max(MT_depth_scaledmappedreads),
            MG_max = max(MG_depth_scaledmappedreads))

# Determine active/inactive by MT/MG ratio and MP presence
CUTOFF_RATIO <- 1
profiles <- profiles %>%
  filter(MG_mean > 0) %>%
  mutate(ratio = MT_mean/MG_mean) %>%
  mutate(Active = ifelse(MP_sum >= 2, 1, ifelse(ratio >= CUTOFF_RATIO, 1, 0)))