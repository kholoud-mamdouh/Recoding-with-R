#performing ordination of bins by expressed KOs at different time-points.
# Load libraries
library(data.table)
library(tidyverse)

# Load tables
samples <- read.table("sample_list.txt")
regelist <- read.table("representative_bins.filtered.tax.tsv")

# Load annotations
anno <- readRDS("Output/Reges_annotations.RDS")
ko_def <- read_tsv("ko_definition.tsv")
anno <- left_join(anno, ko_def)

# Load taxonomy
tax <- fread("tax_reges_filt.amphora2sourmash.filtered.tsv", sep = "\t", header = TRUE)
tax[tax == "low confidence assignment"] <- NA
tax[tax == "unassigned"] <- NA

# Filter for specific KOs
has_dgat_list <- anno %>%
  filter(KO %in% c("K00635", "K18851", "K15406", "K11155", "K11160", "K14456")) %>%
  select(binID) %>%
  unique()

# Load KO expression data
fulltib <- readRDS("Output/MT_MP_profiles_params.rds")

# Gather and join data with annotations
fulltib %>%
  select(binID, data) %>%
  mutate(data = map(data, function(x) gather(x, tp, value, -GeneID))) %>%
  mutate(expr_genes_anno = map(data, function(x) left_join(x, anno, by = "binID")))

# Perform ordination of bins by expressed KOs
''' principal component analysis (PCA), multidimensional scaling (MDS), 
or non-metric multidimensional scaling (NMDS).'''
# Load libraries
library(vegan)

# Load KO expression data
data <- read.csv("KO_expression_data.csv")

# Transform data into matrix format
matrix_data <- as.matrix(data[, -1])

# Standardize data
standardized_data <- scale(matrix_data)

# Perform PCA
pca_result <- prcomp(standardized_data)

# Plot the results
plot(pca_result$x[, 1:2], col = data$time_point, pch = 21)
legend("topleft", legend = unique(data$time_point), col = unique(data$time_point), pch = 21)

# Output is saved in "Output/expressed_profiles_cluster.RDS"
