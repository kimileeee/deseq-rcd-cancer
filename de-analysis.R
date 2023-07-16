## 
library(TCGAbiolinks)
library(tidyverse)

# get a list of projects
gdcprojects <- getGDCprojects()

# Querying Colorectal Cancer Data
query <- GDCquery(
  project = "TCGA-COAD",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  experimental.strategy = "RNA-Seq",
  workflow.type = "STAR - Counts",
  access = open
)
getResults(query)
GDCdownload(query)

# prepare count data
tcga_coad_data <- GDCprepare(query_TCGA, summarizeExperiment = TRUE)
# unstanded, stranded_first, stranded_second, tpm_unstrand, fpkm_unstrand, fpkm_uq_unstrand
coad_matrix <- assay(tcga_coad_data, 'unstranded')

# build a query to get DNA




library(DESeq2)
library(tidyverse)
library(airway)

counts_data <- read_delim( file = "counts_data.tsv", delim = "\t",
                   escape_double = FALSE,
                   trim_ws = TRUE, skip = 1  ) %>%
               slice(5:n())
head(counts_data)

