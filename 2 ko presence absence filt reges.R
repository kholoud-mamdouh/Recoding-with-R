''' generating a presence-absence matrix of KEGG annotations,
and then it performs multidimensional scaling (MDS) 
on the resulting distance matrix to generate an embedding of the data.'''
# Load libraries
library(tidyverse)
library(vegan)

# Load data
regelist <- read_tsv("representative_bins.filtered.tax.tsv", col_names = "ReGeID")
all_genes <- fread("All_genes_CMPs.tsv")

# Filter and select KOs per bin
rege_ko <- all_genes %>% 
  filter(binID %in% regelist$ReGeID) %>%
  select(binID, KO) %>%
  mutate(KO = strsplit(as.character(KO), ";") %>% unnest() %>% filter(!is.na(KO)) %>% gsub( "ko:","", KO))

# Convert to presence-absence matrix
rege_ko.dup$value = 1
pa_mat = spread(rege_ko.dup, KO, value, fill = 0)

# Calculate distance and embedding
d.veg <- vegdist(pa_mat, method="jaccard", binary=T)
fit.mds <- cmdscale(d.veg, eig=T)
embedding_mds.df <- data.frame(x = fit.mds$points[,1], y = fit.mds$points[,2], genome = row.names(pa_mat))
