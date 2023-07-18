## Downloading TCGA Data in R
#. For downloading colorectal and breast cancer gene data from the Genomic Data 
#. Commons (GDC) portal onto local machine
#. Documentation: https://bioconductor.org/packages/release/bioc/html/TCGAbiolinks.html

# Pre-requisite: Installing TCGAbiolinks
if (!requireNamespace("BiocManager", quietly=TRUE))
  install.packages("BiocManager")
BiocManager::install("TCGAbiolinks")

# Libraries
library(TCGAbiolinks)
library(dplyr)
library(DT) 

#. Guide on searching the GDC database: 
#. https://bioconductor.org/packages/release/bioc/vignettes/TCGAbiolinks/inst/doc/query.html

## Colorectal Cancer Query
query <- GDCquery(
  project = "TCGA-COAD",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  experimental.strategy = "RNA-Seq",
  workflow.type = "STAR - Counts",
  access = "open",
  barcode = c("TCGA-D5-6540-01A", "TCGA-AA-3525-11A", "TCGA-AA-3525-01A", "TCGA-AA-3815-01A", "TCGA-D5-6923-01A")
)
getResults(query)
GDCdownload(query)
data.colorectal <- GDCprepare(query)


## Breast Cancer Query
query <- GDCquery(
  project = "TCGA-BRCA",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  experimental.strategy = "RNA-Seq",
  workflow.type = "STAR - Counts",
  access = "open"
)
getResults(query)
GDCdownload(query)
data.breast <- GDCprepare(query)
