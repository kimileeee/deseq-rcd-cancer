## Downloading TCGA Data in R
#. For downloading colorectal and breast cancer gene data from the Genomic Data 
#. Commons (GDC) portal onto local machine


# Pre-requisite: Installing TCGAbiolinks
if (!requireNamespace("BiocManager", quietly=TRUE))
  install.packages("BiocManager")
BiocManager::install("TCGAbiolinks")

# Libraries
library(TCGAbiolinks)
library(dplyr)
library(DT)

# Colorectal Cancer Query
query <- GDCquery(
  project = "TCGA-COAD",
  data.category = "DNA methylation",
  barcode = c("TCGA-06-0122","TCGA-14-1456"),
  platform = "Illumina Human Methylation 27",
  legacy = TRUE
)
GDCdownload(query)
data.colorectal <- GDCprepare(query)

# Breast Cancer Query
query <- GDCquery(
  project = "TCGA-COAD",
  data.category = "DNA methylation",
  barcode = c("TCGA-06-0122","TCGA-14-1456"),
  platform = "Illumina Human Methylation 27",
  legacy = TRUE
)
GDCdownload(query)
data.colorectal <- GDCprepare(query)