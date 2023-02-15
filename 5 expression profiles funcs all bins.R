#analyzing the presence absence of KEGG annotations in all bins.
# Load necessary libraries
library(tidyverse)

# load representative bins and funcs assignment
regelist_rmags <- read_tsv("representative_bins.filtered.tax.tsv", col_names = "ReGeID")
funcs <- readRDS("kmeans_ko_presabs.RDS")[[2]]
colnames(funcs) = c("binID","FunC_ID")
funcs$FunC_ID = as.factor(funcs$FunC_ID)

# load checkm result and taxonomy
checkm_pgolb <- read_tsv("checkM_results_PGOL.tsv") %>% 
  mutate(binID=gsub(".fa$|.fasta|.fna$","", `Bin Id`))
tax.combined <- "tax_reges_all.amphora2sourmash.filtered.tsv"
tax.tab <- fread(tax.combined,sep="\t",header=T)
tax <- apply(tax.tab,2,function(x) gsub("\\(\\d+.+\\d+\\)","",x))
tax <- as.data.frame(tax)
tax[tax=="low confidence assignment"] <- NA
tax[tax=="unassigned"] <- NA

# load KO annotation files
files <- list.files("Databases/Annotations/KEGG", full.names = T)
all_genes <- tibble(files) %>% 
  mutate(KO_tab = map(files, read_tsv))

# extract bin_KO
bin_ko <- all_genes %>%
  unnest(KO_tab) %>%
  mutate(binID = gsub("([^_]+_[^_]+)_PROKKA_\\d+","\\1", Gene)) %>%
  filter(binID %in% checkm_pgolb[checkm_pgolb$Completeness>50,]$binID) %>% 
  select(binID, KO=ID) %>% 
  distinct()
