''' gene expression, depth of coverage in mitochondrial and chloroplast genomes,
taxonomy, sample information, and metabolomics data.'''

# Load necessary libraries
library(data.table)
library(tidyverse)

# Read in configuration file
source("config.R")

# Set output directory for processed data
outdir <- paste(RESULTSDIR,"PopulationAnalysis/Calculations/Analysed",sep="/")

# Set file paths for mitochondrial and chloroplast depth data
mg.out <- paste(outdir,"mg_depth_genelevel_bins.normalized_totalreads.tsv",sep="/")
mt.out <- paste(outdir,"mt_depth_genelevel_bins.normalized_totalreads.tsv",sep="/")

# Read in mitochondrial and chloroplast depth data
mg.depth <- fread(mg.out)
mt.depth <- fread(mt.out)

# Set file path for gene expression data
ko.file <- paste(RESULTSDIR,"reges_filtered.ko_cmp_annotation.deseq_len.tsv",sep="/")

# Read in gene expression data
ko_tab <- fread(ko.file, sep="\t",header=T)

# Set file path for filtered bin data
filt.list <- read.table(paste(REPODIR,"aux_datafiles","representative_bins.filtered.tax.tsv",sep="/"))

# Set file path for taxonomy data
tax.combined <- paste(RESULTSDIR,"tax_reges_filt.amphora2sourmash.filtered.tsv",sep="/")

# Read in taxonomy data, remove scores and low-confidence assignments
tax.tab <- fread(tax.combined,sep="\t",header=T)
tax <- apply(tax.tab,2,function(x) gsub("\\(\\d+.+\\d+\\)","",x))
tax <- as.data.frame(tax)
tax[tax=="low confidence assignment"] <- NA
tax[tax=="unassigned"] <- NA

# Set file path for sample data
samplelist <- paste(REPODIR,"aux_datafiles/sample_list.txt",sep="/")

# Read in sample data and process date information
samples <- read.table(samplelist,sep="\t")
colnames(samples) <- c("sample","Date")
rownames(samples) <- samples$sample
samples <- samples[-c(1,2),]
samples <- samples[-1]
samples$Month <- gsub("\\d+-(\\d+)-\\d+","\\1",samples$Date)
samples$Year <- gsub("(\\d+)-\\d+-\\d+","\\1",samples$Date)
samples$SampleID <- rownames(samples)
samples <- arrange(samples,Date)

# Set file path for physico-chemical data
pc.file <- paste(RESULTSDIR,"/Databases/PhysicoChemical/pc_params.processed_sampling_dates.tsv",sep="/")

# Read in physico-chemical data
pcpar <- fread(pc.file)
pcpar.long <- dcast(pcpar,Date~variable,value.var="value")

# Set file path for metabolomics data
metfile <- paste(
  