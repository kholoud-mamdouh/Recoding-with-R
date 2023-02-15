''' loads the data from two RDS files: "kmeans_ko_presabs.RDS" and 
"Reges_annotations.RDS" that contains clustering results and annotations
of the rMAGs respectively.
Then it selects the columns of interest and renames them in the dataframe.
It filters the annotations to include only the filtered rMAGs that belong to
a FunC using a left join.
It then unnests any multiple entries in the annotations,
groups the data by FunC and KO, and counts the number of occurrences of each KO 
in each FunC.
It calculates several ratios,
such as the ratio of the number of KOs present in a FunC 
to the total number of KOs, and the ratio of the number of bins in which 
a KO is present to the total number of bins'''

# Load libraries
library(data.table)
library(stringr)
library(tidyverse)
library(vegan)
library(broom)
library(here)
library(readxl)

# Load data
clustres <- readRDS("kmeans_ko_presabs.RDS")

# Select columns and rename
func_ids <- clustres[[2]] %>% 
  rename(binID = "group", FunC_ID = "clusterID") %>% 
  mutate(FunC_ID = as.factor(FunC_ID))

# Load annotations
anno <- readRDS("Reges_annotations.RDS")

# Filter annotations by func_ids
anno <- left_join(anno, func_ids)

# Unnest multiple entries and group by FunC and KO
anno_unnest <- anno %>%
  filter(!is.na(FunC_ID)) %>%
  mutate(KO = gsub("ko:","",KO)) %>%
  mutate(KO = strsplit(KO,";")) %>%
  unnest(KO)

# Group by FunC and KO and count
ko_counts_total <- table(anno_unnest$KO) %>% as.data.frame()
colnames(ko_counts_total) = c("KO","Freq_KOs")

funcs_num <- table(func_ids$FunC_ID) %>% as.data.frame()
colnames(funcs_num) = c("FunC_ID","Freq_bins")

binfuncko <- anno_unnest %>% 
  group_by(FunC_ID, KO, binID) %>%
  summarise(count=length(KO))

# Join dataframes to calculate ratios
number_funcKO <- binfuncko %>% 
  filter(!is.na(KO)) %>%
  group_by(KO, FunC_ID) %>%
  summarise(count_perfunc = sum(count), presence_bins = length(unique(binID)))

final <- number_funcKO %>%
  left_join(ko_counts_total) %>%
  left_join(funcs_num) %>%
  mutate(ratio_KO_presence = count_perfunc / Freq_KOs) %>% 
  mutate(occurence_bins = presence_bins / Freq_bins)
